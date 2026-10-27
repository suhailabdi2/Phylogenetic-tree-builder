package ncbi

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"
	"log/slog"
	"net/http"
	"net/url"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"sync"
	"time"

	"github.com/suhailabdi2/Phylogenetic-tree-builder/fasta"
)

const (
	defaultBaseURL     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
	defaultMinInterval = 350 * time.Millisecond // a little over 1/3s => under 3 req/s
	defaultTimeout     = 30 * time.Second
	maxBodyBytes       = 256 << 20 // refuse responses larger than 256 MiB
	toolName           = "phylo"
)

var (
	// ErrInvalidAccession means the accession failed local validation; no
	// request was made.
	ErrInvalidAccession = errors.New("ncbi: invalid accession")

	// ErrNotFound means NCBI answered successfully but returned nothing.
	ErrNotFound = errors.New("ncbi: accession not found (empty response)")

	// ErrBadResponse means NCBI returned something that isn't the single
	// FASTA record we asked for (an error message, the wrong record, ...).
	ErrBadResponse = errors.New("ncbi: unexpected response")
)

// HTTPError is returned when NCBI answers with a non-200 status.
type HTTPError struct {
	Status  int
	Snippet string // start of the response body, for debugging
}

func (e *HTTPError) Error() string {
	msg := fmt.Sprintf("ncbi: http status %d", e.Status)
	if e.Status == http.StatusTooManyRequests {
		msg += " (rate limited)"
	}
	if e.Snippet != "" {
		msg += fmt.Sprintf(": %q", e.Snippet)
	}
	return msg
}

// accessionRE matches the accession shapes we expect:
//
//	NC_012920.1      RefSeq: 2 letters, '_', digits
//	NZ_CP012345.1    RefSeq WGS: 2 letters, '_', letters, digits
//	CM018432.1       GenBank: letters then digits
//	VOAJ01024641.1   GenBank WGS: 4 letters then digits
//
// An optional ".version" suffix is allowed. The pattern is deliberately
// strict: the accession ends up in a URL and in a file name, so nothing like
// '/', '&' or '..' may get through.
var accessionRE = regexp.MustCompile(`^(?:[A-Z]{2}_[A-Z]{0,4}|[A-Z]{1,6})[0-9]{5,}(?:\.[0-9]+)?$`)

func normalizeAccession(s string) (string, error) {
	acc := strings.ToUpper(strings.TrimSpace(s))
	if !accessionRE.MatchString(acc) {
		return "", fmt.Errorf("%w: %q", ErrInvalidAccession, s)
	}
	return acc, nil
}

// Client fetches sequences from NCBI. It is safe for concurrent use.
type Client struct {
	cacheDir    string
	email       string
	baseURL     string
	httpClient  *http.Client
	minInterval time.Duration

	mu   sync.Mutex // guards next
	next time.Time  // earliest time the next request may be sent
}

// Option customises a Client.
type Option func(*Client)

// WithBaseURL overrides the efetch endpoint (used by tests).
func WithBaseURL(u string) Option { return func(c *Client) { c.baseURL = u } }

// WithHTTPClient replaces the default HTTP client (30s timeout).
func WithHTTPClient(h *http.Client) Option { return func(c *Client) { c.httpClient = h } }

// WithMinInterval sets the minimum gap between network requests.
// Use a smaller value only if you have an NCBI API key.
func WithMinInterval(d time.Duration) Option { return func(c *Client) { c.minInterval = d } }

// NewClient returns a Client that caches downloads in cacheDir (created on
// first write). email is sent to NCBI so they can contact you if your usage
// causes problems; it may be empty.
func NewClient(cacheDir, email string, opts ...Option) *Client {
	c := &Client{
		cacheDir:    cacheDir,
		email:       email,
		baseURL:     defaultBaseURL,
		httpClient:  &http.Client{Timeout: defaultTimeout},
		minInterval: defaultMinInterval,
	}
	for _, o := range opts {
		o(c)
	}
	return c
}

// Fetch returns the sequence for one accession, using the disk cache when it
// can. The accession may be given with or without a version
// ("NC_012920" or "NC_012920.1"); the cache is keyed by what you pass.
func (c *Client) Fetch(ctx context.Context, accession string) (fasta.Record, error) {
	acc, err := normalizeAccession(accession)
	if err != nil {
		return fasta.Record{}, err
	}

	if rec, ok := c.readCache(acc); ok {
		return rec, nil
	}

	body, err := c.download(ctx, acc)
	if err != nil {
		return fasta.Record{}, fmt.Errorf("fetch %s: %w", acc, err)
	}
	rec, err := parseSingle(acc, body)
	if err != nil {
		return fasta.Record{}, fmt.Errorf("fetch %s: %w", acc, err)
	}

	// A failed cache write shouldn't lose a sequence we already paid for,
	// but it shouldn't be silent either.
	if err := c.writeCache(acc, body); err != nil {
		slog.Warn("ncbi: could not write cache", "accession", acc, "err", err)
	}
	return rec, nil
}

// FetchAll fetches each accession in order and returns the records in the
// same order. It stops at the first failure. Requests run one at a time,
// which is what NCBI's rate limit wants anyway.
func (c *Client) FetchAll(ctx context.Context, accessions []string) ([]fasta.Record, error) {
	// Validate everything first so a typo in the last accession doesn't
	// cost us the downloads before it.
	for _, acc := range accessions {
		if _, err := normalizeAccession(acc); err != nil {
			return nil, err
		}
	}

	out := make([]fasta.Record, 0, len(accessions))
	for _, acc := range accessions {
		rec, err := c.Fetch(ctx, acc)
		if err != nil {
			return nil, err
		}
		out = append(out, rec)
	}
	return out, nil
}

// --- network ---------------------------------------------------------------

func (c *Client) download(ctx context.Context, acc string) ([]byte, error) {
	if err := c.wait(ctx); err != nil {
		return nil, err
	}

	q := url.Values{}
	q.Set("db", "nuccore")
	q.Set("id", acc)
	q.Set("rettype", "fasta")
	q.Set("retmode", "text")
	q.Set("tool", toolName)
	if c.email != "" {
		q.Set("email", c.email)
	}

	req, err := http.NewRequestWithContext(ctx, http.MethodGet, c.baseURL+"?"+q.Encode(), nil)
	if err != nil {
		return nil, err
	}
	req.Header.Set("User-Agent", toolName)

	resp, err := c.httpClient.Do(req)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	// Read one byte past the limit so we can tell "exactly at" from "over".
	body, err := io.ReadAll(io.LimitReader(resp.Body, maxBodyBytes+1))
	if err != nil {
		return nil, err
	}
	if len(body) > maxBodyBytes {
		return nil, fmt.Errorf("response larger than %d bytes", maxBodyBytes)
	}

	if resp.StatusCode != http.StatusOK {
		return nil, &HTTPError{Status: resp.StatusCode, Snippet: snippet(body)}
	}
	return body, nil
}

// wait blocks until this request may be sent, or ctx is done. Each caller
// reserves the next free slot, so concurrent callers queue up in order.
func (c *Client) wait(ctx context.Context) error {
	c.mu.Lock()
	now := time.Now()
	slot := c.next
	if slot.Before(now) {
		slot = now
	}
	c.next = slot.Add(c.minInterval)
	c.mu.Unlock()

	d := time.Until(slot)
	if d <= 0 {
		return ctx.Err()
	}
	t := time.NewTimer(d)
	defer t.Stop()
	select {
	case <-t.C:
		return nil
	case <-ctx.Done():
		return ctx.Err()
	}
}

// --- parsing ---------------------------------------------------------------

// parseSingle parses data and checks it is exactly the record we asked for.
func parseSingle(acc string, data []byte) (fasta.Record, error) {
	if len(bytes.TrimSpace(data)) == 0 {
		return fasta.Record{}, ErrNotFound
	}
	recs, err := fasta.Parse(bytes.NewReader(data))
	if err != nil {
		return fasta.Record{}, fmt.Errorf("%w: %v (response began %q)", ErrBadResponse, err, snippet(data))
	}
	if len(recs) != 1 {
		return fasta.Record{}, fmt.Errorf("%w: got %d records, want 1", ErrBadResponse, len(recs))
	}
	rec := recs[0]
	// NCBI may add a version to an unversioned request ("NC_012920" ->
	// "NC_012920.1"), but it must be the same accession. This also stops a
	// wrong record from being written to the cache under this name.
	if rec.ID != acc && !strings.HasPrefix(rec.ID, acc+".") {
		return fasta.Record{}, fmt.Errorf("%w: asked for %s, got %s", ErrBadResponse, acc, rec.ID)
	}
	return rec, nil
}

func snippet(b []byte) string {
	const n = 120
	s := strings.TrimSpace(string(b))
	if len(s) > n {
		s = s[:n] + "..."
	}
	return s
}

// --- cache -----------------------------------------------------------------

func (c *Client) cachePath(acc string) string {
	// acc has passed accessionRE, so it is safe as a file name.
	return filepath.Join(c.cacheDir, acc+".fasta")
}

// readCache returns the cached record, or ok=false for a miss. A corrupt or
// mismatched cache file counts as a miss and will be overwritten.
func (c *Client) readCache(acc string) (rec fasta.Record, ok bool) {
	data, err := os.ReadFile(c.cachePath(acc))
	if err != nil {
		if !errors.Is(err, os.ErrNotExist) {
			slog.Warn("ncbi: could not read cache", "accession", acc, "err", err)
		}
		return fasta.Record{}, false
	}
	rec, err = parseSingle(acc, data)
	if err != nil {
		slog.Warn("ncbi: ignoring bad cache file", "accession", acc, "err", err)
		return fasta.Record{}, false
	}
	return rec, true
}

// writeCache writes atomically (temp file + rename) so a crash mid-write
// can't leave a half-written file that later looks like a valid cache entry.
func (c *Client) writeCache(acc string, data []byte) error {
	if err := os.MkdirAll(c.cacheDir, 0o755); err != nil {
		return err
	}
	tmp, err := os.CreateTemp(c.cacheDir, acc+".*.tmp")
	if err != nil {
		return err
	}
	tmpName := tmp.Name()
	if _, err := tmp.Write(data); err != nil {
		tmp.Close()
		os.Remove(tmpName)
		return err
	}
	if err := tmp.Close(); err != nil {
		os.Remove(tmpName)
		return err
	}
	if err := os.Rename(tmpName, c.cachePath(acc)); err != nil {
		os.Remove(tmpName)
		return err
	}
	return nil
}
