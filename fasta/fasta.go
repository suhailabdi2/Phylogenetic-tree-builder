// Package fasta parses nucleotide sequences in FASTA format.
//
// A FASTA file is a list of records. Each record starts with a header line
// beginning with '>' followed by one or more lines of sequence:
//
//	>NC_012920.1 Homo sapiens mitochondrion, complete genome
//	GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT
//	TTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTG
package fasta

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
)

// Record is a single named sequence from a FASTA file.
type Record struct {
	ID          string // first word of the header, e.g. "NC_012920.1"
	Description string // rest of the header line (may be empty)
	Seq         string // upper-case sequence with line breaks removed
}

// ParseError reports a problem at a specific line of the input.
type ParseError struct {
	Line int
	Msg  string
}

func (e *ParseError) Error() string {
	return fmt.Sprintf("fasta: line %d: %s", e.Line, e.Msg)
}

// ErrNoRecords is returned when the input contains no sequences at all.
var ErrNoRecords = errors.New("fasta: no records found")

// validResidues are the characters allowed in a sequence: the IUPAC
// nucleotide codes plus '-' for alignment gaps. Input is upper-cased first.
const validResidues = "ACGTUNRYKMSWBDHV-"

var validTable [256]bool

func init() {
	for i := 0; i < len(validResidues); i++ {
		validTable[validResidues[i]] = true
	}
}

// Parse reads every record from r.
//
// It rejects: sequence data before the first header, empty headers, empty
// records, duplicate IDs, and any character that is not a valid nucleotide
// code. Blank lines and ';' comment lines are ignored, and both LF and CRLF
// line endings are accepted.
func Parse(r io.Reader) ([]Record, error) {
	br := bufio.NewReader(r)

	var (
		records []Record
		cur     *Record
		curLine int // line number of cur's header
		sb      strings.Builder
		seenAt  = make(map[string]int) // ID -> header line, for duplicate detection
		lineNo  int
	)

	// finish closes the record currently being built.
	finish := func() error {
		if cur == nil {
			return nil
		}
		if sb.Len() == 0 {
			return &ParseError{Line: curLine, Msg: fmt.Sprintf("record %q has no sequence", cur.ID)}
		}
		cur.Seq = sb.String()
		records = append(records, *cur)
		cur = nil
		sb.Reset()
		return nil
	}

	for {
		// ReadString (unlike bufio.Scanner) has no line-length limit, which
		// matters because genomes are often stored on a single huge line.
		raw, readErr := br.ReadString('\n')
		if readErr != nil && readErr != io.EOF {
			return nil, fmt.Errorf("fasta: read: %w", readErr)
		}

		if raw != "" {
			lineNo++
			line := strings.TrimRight(raw, "\r\n")

			switch {
			case strings.HasPrefix(line, ";"):
				// old-style comment line: ignore

			case strings.HasPrefix(line, ">"):
				if err := finish(); err != nil {
					return nil, err
				}
				header := strings.TrimSpace(line[1:])
				if header == "" {
					return nil, &ParseError{Line: lineNo, Msg: "empty header after '>'"}
				}
				id, desc := header, ""
				if i := strings.IndexAny(header, " \t"); i >= 0 {
					id, desc = header[:i], strings.TrimSpace(header[i+1:])
				}
				if first, dup := seenAt[id]; dup {
					return nil, &ParseError{
						Line: lineNo,
						Msg:  fmt.Sprintf("duplicate ID %q (first seen on line %d)", id, first),
					}
				}
				seenAt[id] = lineNo
				cur = &Record{ID: id, Description: desc}
				curLine = lineNo

			default:
				seq := strings.TrimSpace(line)
				if seq == "" {
					break // blank line
				}
				if cur == nil {
					return nil, &ParseError{Line: lineNo, Msg: "sequence data before first '>' header"}
				}
				for i := 0; i < len(seq); i++ {
					c := seq[i]
					if c >= 'a' && c <= 'z' {
						c -= 'a' - 'A'
					}
					if !validTable[c] {
						return nil, &ParseError{
							Line: lineNo,
							Msg:  fmt.Sprintf("invalid character %q at column %d in record %q", seq[i], i+1, cur.ID),
						}
					}
					sb.WriteByte(c)
				}
			}
		}

		if readErr == io.EOF {
			break
		}
	}

	if err := finish(); err != nil {
		return nil, err
	}
	if len(records) == 0 {
		return nil, ErrNoRecords
	}
	return records, nil
}

// ParseFile opens the named file and parses it with Parse.
func ParseFile(path string) ([]Record, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	records, err := Parse(f)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return records, nil
}
