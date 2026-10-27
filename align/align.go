package align

import (
	"errors"
	"strings"
)

// ErrGapInInput is returned when either input already contains the gap marker.
var ErrGapInInput = errors.New("input sequence contains a gap")

// Scoring assigns scores to matches, mismatches, and gaps.
type Scoring struct {
	Match, Mismatch, Gap int
}

// DefaultScoring uses +1 for a match, -1 for a mismatch, and -2 for a gap.
var DefaultScoring = Scoring{Match: 1, Mismatch: -1, Gap: -2}

// Alignment contains two equal-length aligned sequences and their score.
type Alignment struct {
	A, B  string
	Score int
}

// NeedlemanWunsch computes an optimal global alignment using linear gap costs.
func NeedlemanWunsch(a, b string, scoring Scoring) (Alignment, error) {
	if strings.ContainsRune(a, '-') || strings.ContainsRune(b, '-') {
		return Alignment{}, ErrGapInInput
	}

	cols := len(b) + 1
	trace := make([]byte, (len(a)+1)*cols)
	previous := make([]int, cols)
	current := make([]int, cols)
	for j := 1; j <= len(b); j++ {
		previous[j] = j * scoring.Gap
		trace[j] = 'l'
	}

	for i := 1; i <= len(a); i++ {
		current[0] = i * scoring.Gap
		trace[i*cols] = 'u'
		for j := 1; j <= len(b); j++ {
			idx := i*cols + j
			match := scoring.Mismatch
			if a[i-1] == b[j-1] {
				match = scoring.Match
			}
			best, direction := previous[j-1]+match, byte('d')
			if candidate := previous[j] + scoring.Gap; candidate > best {
				best, direction = candidate, 'u'
			}
			if candidate := current[j-1] + scoring.Gap; candidate > best {
				best, direction = candidate, 'l'
			}
			current[j], trace[idx] = best, direction
		}
		previous, current = current, previous
	}

	alignedA := make([]byte, 0, len(a)+len(b))
	alignedB := make([]byte, 0, len(a)+len(b))
	for i, j := len(a), len(b); i > 0 || j > 0; {
		switch trace[i*cols+j] {
		case 'd':
			alignedA = append(alignedA, a[i-1])
			alignedB = append(alignedB, b[j-1])
			i, j = i-1, j-1
		case 'u':
			alignedA = append(alignedA, a[i-1])
			alignedB = append(alignedB, '-')
			i--
		case 'l':
			alignedA = append(alignedA, '-')
			alignedB = append(alignedB, b[j-1])
			j--
		}
	}
	for i, j := 0, len(alignedA)-1; i < j; i, j = i+1, j-1 {
		alignedA[i], alignedA[j] = alignedA[j], alignedA[i]
		alignedB[i], alignedB[j] = alignedB[j], alignedB[i]
	}

	return Alignment{
		A:     string(alignedA),
		B:     string(alignedB),
		Score: previous[len(b)],
	}, nil
}
