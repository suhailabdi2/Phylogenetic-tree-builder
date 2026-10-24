package main

import (
	"fmt"
	"os"

	"github.com/suhailabdi2/Phylogenetic-tree-builder/fasta"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Fprintln(os.Stderr, "usage: fastainfo <file.fasta> [more files...]")
		os.Exit(2)
	}
	for _, path := range os.Args[1:] {
		records, err := fasta.ParseFile(path)
		if err != nil {
			fmt.Printf("FAIL  %v\n", err)
			continue
		}
		fmt.Printf("OK    %s: %d record(s)\n", path, len(records))
		for _, r := range records {
			fmt.Printf("      %-14s %8d bp  %s\n", r.ID, len(r.Seq), r.Description)
		}
	}
}
