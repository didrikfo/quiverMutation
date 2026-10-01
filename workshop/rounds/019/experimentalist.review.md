# Review of workshop/rounds/019/experimentalist.md

referee: theorist · round: 019
verdict: minor revision

## Reproduction

Re-ran `--e078` (3 s): output matches exactly (75 tilting / 6 non-tilting; total dim ker histogram 3/2/1; per-i 22 zero, 10 one; 0 disagreements). Did not re-run n = 5 (82 s) or the two 8-minute capped runs; read the committed outputs instead. Their tails are internally consistent with the table: 30 300 + 0, 83 591 + 907, 53 502 + 143 steps; per-i counts satisfy zero = (n-1)*non-tilting - ones (3628 = 4*907, 715 = 5*143 minus ones, i.e. consistent); the sums 75+30 300+83 591+53 502 = 167 468 and 6+907+143 = 1 056 agree with the claim. The n = 6, 7 numbers are time-capped, so a re-run would not reproduce them exactly (author says so).

## True?

The numbers are as stated. Two problems of statement, not of data.

1. Wording "Cartan discrepancy = -dim ker ... on every gate-admitted step tested" is fine, but the claim "the converse (tilting => dim ker 0) is the new content" is not new content. E-093 already records that the discrepancy is row k off-diagonal with entry (k,i) = -dim ker g_i, and that congruence fails exactly where tiltingPlus fails (no disagreements). A tilting step is congruent, so the discrepancy is 0, so dim ker = 0 by the E-093 identity. The 0-of-167 000 "tilting with dim ker > 0" is that same identity read on tilting steps; it is a consistency check of E-093's code path (same script plus tally), not an independent test.
2. "dim ker g_i never above 1" holds only for the 1 050 walk steps plus the E-078 family (where g_i is also <= 1; totals 2 and 3 come from several i, each 1). The title says "on the guarded walks", correct, but "each rejecting parent has exactly one non-tilting vertex" is about vertices of the parent, and the sum over i != k at one vertex gives total 1 only on the walks. These are different quantities and the Claim runs them together once ("single entry -1" vs "dim ker never above 1"). Minor.

No counterexample found. The author's own caveat (parents all A5-shaped per E-084) is the real limit: the "<= 1" observation is a statement about one parent shape, not about the gate.

## New?

grep of `research/` for "dim ker": E-093 (EXPERIMENTS.md line 9/10), whose Limits say the dim ker distribution and distinct-parent count were not reported. So the histogram and the distinct-parent count (907 and 143 distinct rejecting parents, one non-tilting vertex each) are new; the tilting <=> dim ker 0 equivalence is E-093 restated. H-015 and E-084 related as cited.

## Evidenced?

Mostly. Counts, caps, command lines and output files are given. Missing: (a) the n = 6 versus E-093 difference (907 vs 696) is attributed to machine load without a check; the same cap in steps would settle it; (b) no per-parent shape check, which the author admits and which is exactly what decides whether "dim ker <= 1" is about the gate or about A5; (c) n = 6 classes 1-3, n = 7 classes 1+ not run, stated.

## Required for acceptance

1. Remove "the converse is the new content" or say plainly it follows from E-093's identity on its data; keep the histogram and distinct-parent counts as the new content.
2. Either do the A5-shape check of the 1 050 rejecting parents (the author's own "Next") or state in the Claim that the "<= 1" observation is conditional on the walks' parents being A5-shaped, untested here.
3. Separate "number of non-tilting vertices per parent (always 1)" from "dim ker per vertex (always 1 on the walks)" in the Claim.
