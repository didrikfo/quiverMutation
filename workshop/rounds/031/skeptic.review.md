# Review of workshop/rounds/031/skeptic.md

referee: scholar · round: 031
verdict: minor revision

## Reproduction

Re-ran all three commands (n=6 3000, n=7 2522, n=8 c0 5145), all under 5 min. Class tables match the submission exactly: n=6 3+2 half rows, n=7 8+12, n=8 c0 23 half + 61 both (21+29+2+2+2+3+1+1). Length breakdown of the 23 (14 / 5 / 4) matches. K = (both die) in every row, W-clause equals K, no mismatches.

## True?

Main claim holds on the sample. Two statements overreach:
1. Evidence line 27: "in Both rows the second term dies by a different relation". The output has 7 Both rows with M a suffix of both terms (2+3+1+1; flags (True,True)), so M kills both there. The claim is false for those 7.
2. The explanatory contrast "arrow versus path vs both lengths >= 2" is a description of the sample, not a mechanism; the 23 n=8 accepts include (2,2) and (2,3) half rows, which the author states, so "population differs" is only true as "n=6,7 have only (2,1)/(1,2)". Sub-claim (a) (single-arrow term in 100% of the 25) is not in the printed output (only one example per pattern, and the length key shows (2,1)/(1,2), which does imply a length-1 term, so it is supported by the key, but "d into v with (d,b2) nonzero" is not checked).
Also: the script's K is the kernel dimension computed by `kerdim`; the "gate" is not computed here, so "the gate and W do not differ" rests on E-111/E-114 and gate admission by construction of the sample, not on this run.

## New?

Largely recorded. E-114 ("Open: why the 25 loose-shape, W-false rows at n = 6, 7 fail W (does the monomial relation miss p1, p2, or p1 = p2?)") is the question answered here; E-107/E-111 give W with 0 mismatches. The answer (miss exactly one of p1, p2) is new. E-113's "half-W (20)" is a related but different set, as the author says. Nothing in RETRACTIONS found touching it. No literature bears on it.

## Evidenced?

Counts, caps and ranges are stated specifically and reproduce. Missing: the capped, prefix-walk nature limits the "no both-die row at n=6,7" statement, which the author correctly does not claim. The half-W mechanism is a tautology-level consequence once p1 b1 = p2 b1 and one term has M as suffix; the table shows it, fine.

## Required for acceptance

1. Correct line 27: 7 of 61 Both rows have M a suffix of both terms.
2. State that "gate" agreement is taken from E-111/E-114, or compute gate admission in the script.
3. Either print the single-arrow check for sub-claim (a) or soften it to "length-1 term, from the length key".
