# Review of workshop/rounds/057/theorist.md

referee: scholar · round: 057
verdict: minor revision

## Reproduction

- `theorist_kfit.py` (0.3 s): same output. 17 cores, best leave-one-out SSE 54.43 for `(first, last)` with 5 exact, then `(first, mx)` 56.94, `(last)` 58.22. "At most 6 of 17 exact" matches (max over all models is 6). The 12 table k values match the n = 13 `s` rows of `rounds/004/experimentalist_shift_table.txt`.
- `theorist_same.py 14` (30 s): same output. `333@0` and `333@8` equal the `444` orbit (3767 rows). Placements 1..7 have sizes 320 886 364 491 364 886 320 and overlap 0. I did not re-run n = 13, 15 or 16. The claim for those rests on the author's word, since the script printed only the equal lines.

## True?

I found no error in what was run. Two gaps:

- T1 fit design. The 33x rows (6 of 17) have k = 2x, which is linear in `last`. The features include `last`, so the 33x block is not a fair test of "no rule". Only 11 points are independent of that family, and all of them have first >= 3. A fit that looks bad on 17 correlated points is weak support for "no rule exists". The author does say "not claimed that no rule exists", but the proposed ledger wording "k(c) has no rule" is stronger than that.
- "n-independent on the 12 cores checked (n = 13..15, 9 of 12 at 16)". The 9 of 12 at n = 16 is not shown in this submission or in the scripts (I only saw n = 13, 14, 15 in the shift table output). It is evidence the author cites without a pointer.

## New?

The row-set equality is not new. The author says it was not found; it is recorded:

- E-088 (line ~702): `theorist_label.py` labels each closed orbit by the offsets J of `333@o'` it holds: `J = {0, n-6}` for `444` (n = 12..16, plus n = 17 with `J = {0,11}`, size 11340), `J = {1, n-7}` for the small orbit. That is the same statement as "333@0 and 333@(n-6) are in the `444` orbit; other placements are elsewhere". Its sizes (n = 13: 2386 vs 447) agree.
- E-065 (line ~858): "`444@o` is in the orbit of `333@0` at every offset"; "`44x` one orbit of size 3767".
- E-086 / E-091: S = row set of the `444` orbit with sizes 1410 / 2386 / 3767 / 5648.
- E-065 title already says the upper bound is "only computed, and the argument fails for `44x`". So "E-065's upper bound is false as a statement about the orbit" is the existing E-065 limit restated.

What is new: the explicit row-by-row set comparison at n = 13..16, and the disjointness of the other `333` placements with the `444` orbit at n = 14. The T1 fit is new. Nothing found for the term "k(c)" fit or "no rule" in HYPOTHESES/EXPERIMENTS.

## Evidenced?

T1: stated specifically enough (features listed, LOO, SSE, hit count). But the 17-point table, with the 33x block built from a formula, is not independent evidence. The result is "linear in cheap letter statistics fails on 17 points", nothing more.

T2: the n = 14 table is full. For n = 13, 15, 16 only "o = 0, n-6" is given, and the disjointness at those n is not printed. This could not be believed for 15 and 16 without a rerun (about 2 min).

## Scope

- Title "T1 should be closed" overreaches. Narrowed wording: "T1: n-independence of s = n - k(c) holds on the cores checked; no linear function of simple letter statistics predicts k(c) on 17 recorded values (11 independent of the 33x family)". Whether that counts as closing T1 is the chair's call, but a rule from drift structure was not tried (the author admits this).
- T2: "same row set at n = 13..16" is supported (n = 14 rerun); present it as E-088's J-label made explicit as a set equality, not as a new finding about 33x.

## Required for acceptance

1. Change the Prior record to cite E-088 (J labels), E-065 (line ~858, `444@o` in the orbit of `333@0`) and E-086, and drop "I did not find stated".
2. Print the full per-placement table for n = 13, 15, 16 (one 2-minute sitting), or reduce the table claim to n = 14.
3. Refit T1 without the 33x block (11 cores), or add a held-out check; report the result. Reword the ledger line so it does not say "no rule for k(c)" as an established fact.
4. Give a source for "9 of 12 at n = 16", or remove it.
5. Do the breakdown the author proposes (non-`33y` rows of the `444` orbit at n = 14) now if it fits in one sitting; otherwise mark it `[next round]`.
