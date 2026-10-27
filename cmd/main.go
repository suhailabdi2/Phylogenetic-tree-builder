// Command fetchdemo downloads sequences from NCBI (with a local cache) and
// prints a one-line summary of each.
//
//	go run ./cmd/fetchdemo -email you@example.com NC_012920 NC_001643 NC_005089
//
// Run it twice: the second run should be instant because everything comes
// from the cache.
package main

import (
	"context"
	"flag"
	"fmt"
	"os"
	"time"

	"github.com/suhailabdi2/Phylogenetic-tree-builder/ncbi"
)

func main() {
	cacheDir := flag.String("cache", ".cache/ncbi", "directory for cached FASTA files")
	email := flag.String("email", "", "contact email sent to NCBI (recommended)")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "usage: fetchdemo [-cache DIR] [-email ADDR] ACCESSION...\n")
		flag.PrintDefaults()
	}
	flag.Parse()
	if flag.NArg() == 0 {
		flag.Usage()
		os.Exit(2)
	}

	client := ncbi.NewClient(*cacheDir, *email)
	ctx, cancel := context.WithTimeout(context.Background(), 2*time.Minute)
	defer cancel()

	start := time.Now()
	recs, err := client.FetchAll(ctx, flag.Args())
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
	for _, r := range recs {
		fmt.Printf("%-14s %8d bp  %s\n", r.ID, len(r.Seq), r.Description)
	}
	fmt.Printf("\n%d sequence(s) in %v\n", len(recs), time.Since(start).Round(time.Millisecond))
}
