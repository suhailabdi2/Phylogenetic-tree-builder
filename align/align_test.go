package align

import (
	"errors"
	"strings"
	"testing"
)

func TestNeedlemanWunsch(t *testing.T) {
	tests := []struct {
		name, a, b, wantA, wantB string
		wantScore                int
	}{
		{"identical", "ACGT", "ACGT", "ACGT", "ACGT", 4},
		{"mismatch", "ACGT", "AGGT", "ACGT", "AGGT", 2},
		{"deletion", "ACGT", "AGT", "ACGT", "A-GT", 1},
		{"one character", "A", "C", "A", "C", -1},
		{"empty first", "", "ACG", "---", "ACG", -6},
		{"empty second", "ACG", "", "ACG", "---", -6},
		{"both empty", "", "", "", "", 0},
		{"no matches", "ACGT", "TGCA", "ACGT", "TGCA", -4},
		{"tie", "AAAC", "AAC", "AAAC", "-AAC", 1},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, err := NeedlemanWunsch(tc.a, tc.b, DefaultScoring)
			if err != nil {
				t.Fatal(err)
			}
			if got.A != tc.wantA || got.B != tc.wantB || got.Score != tc.wantScore {
				t.Errorf("got (%q, %q) score %d; want (%q, %q) score %d",
					got.A, got.B, got.Score, tc.wantA, tc.wantB, tc.wantScore)
			}
		})
	}
}

func TestNeedlemanWunsch_RejectsGaps(t *testing.T) {
	for _, input := range [][2]string{{"AC-GT", "ACGT"}, {"ACGT", "A-GT"}, {"-", "-"}} {
		if _, err := NeedlemanWunsch(input[0], input[1], DefaultScoring); !errors.Is(err, ErrGapInInput) {
			t.Errorf("NeedlemanWunsch(%q, %q) error = %v, want ErrGapInInput", input[0], input[1], err)
		}
	}
}

func TestNeedlemanWunsch_CustomGapScore(t *testing.T) {
	got, err := NeedlemanWunsch("GATTACA", "GCATGCG", Scoring{Match: 1, Mismatch: -1, Gap: -1})
	if err != nil {
		t.Fatal(err)
	}
	if got.Score != 0 || !strings.Contains(got.A+got.B, "-") {
		t.Errorf("got alignment %q / %q, score %d; want score 0 with a gap", got.A, got.B, got.Score)
	}
}
