# Review of workshop/rounds/086/theorist.md

referee: skeptic · round: 086
verdict: accept

## Reproduction

Re-ran the committed script in its final form: `.venv/bin/python -u workshop/rounds/086/theorist_lemmaA.py 60 6 1`, 39 s. Output is byte-identical to theorist_lemmaA_out.txt (diff empty): 300 instances, 13 748 cells, 0 mismatches, 0 overlapping-site words; control 215 instances, 164 with mismatch, 164 with overlapping sites. The script contains no MISMATCH print (grep: only the three summary prints), so the "3 061 filtered lines" are not in either file. Nothing is hidden: the claim table matches the final script.

## True?

Checked the proof, not only the data.
- (1) Overlap of two sites needs a letter in Y (or a theta) that is also in X1. Y and X1 are disjoint and theta is not an arrow, so sites are disjoint. Correct.
- (3) New site after D: needs (z, x) with z in X1, x in Y, but x is in X1; or (y, w) with y in X1, but y is in Y. Both impossible. Also a theta neighbour cannot form a site. Correct. Skeleton preservation gives a direct sum of subcomplexes.
- (4) The piece P (x) S (x) Q is contractible. P has no sites and D=0 on its arrows, so D acts on Q only. The homotopy h (replace xy by theta at the first site) is well defined because the first site is unchanged by D, which is exactly (3). Signs are Koszul and the composition is over R = k^n, with S free on each (x, y) component. No field or characteristic input. Correct.
- (5) The s = 0 skeleta are the normal words with D = 0. Correct.
No gap found. Edge cases I looked for: loops in X1 or Y (fine, only the shared-arrow case overlaps), parallel arrows, a site adjacent to another site's theta, the same arrow used as x in several pairs (distinct thetas, same skeleton pair, harmless).

Minor: the proof is stated for weight-graded finite pieces; T<theta> is infinite-dimensional in total, so "direct sum over skeletons" is of weight-graded pieces. Fine, but should be said.

## New?

Grep of research/ for "Lemma A", Anick, skeleton, tiling, disjoint: only E-235, E-237 (cites it, theorem number UNVERIFIED) and unrelated literature notes. The result is a restatement of the standard no-overlap fact. It is new to the record as a self-contained proof, and E-237's open item is answered by it.

## Evidenced?

Yes. Instance counts, the range (n <= 5, at most 6 arrows, degrees 0..2, weight <= 6, GF(1000003) exact ranks) and the control are all stated. The proof carries the claim, the machine check only corroborates it. The author honestly notes that n = 5 has only 30 composable pairs. Weight <= 6 is irrelevant to a proof valid at all weights.

## Scope

Title and claim match what was checked and what was proved. The claim correctly declines to reach E-235's remaining steps (non-monomial loops at i, d^2 = 0, the spectator step D = 0 after gr_G, journal version). The "Use in E-235" paragraph asserts that step 4 lands in the hypotheses of this lemma. That is an application, not part of this claim, and the author flags it as still to be checked. It should not be recorded as discharged until that check is done.

## Required for acceptance

None required. Suggested for the record, not blocking:
1. When E-235 is annotated, say "caveat discharged for the monomial step, conditional on the step-4 hypotheses check" and not unconditionally.
2. Add one sentence to the proof saying the splitting is per weight.
