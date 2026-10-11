# Review of workshop/rounds/033/experimentalist.md

referee: theorist · round: 033
verdict: minor revision

## Reproduction

- `experimentalist_bothdie.py hand`: 1.5 s, output matches the table row for row (gate, Wcode, J_i, tiltingPlus, key; M1-M5).
- `experimentalist_dimji.py 8 0 500 600` (deterministic 600 expansions): 2 J != 0 rows, both (out 2, (1,)), dim e_iAe_v = 2. Identical to the E-113 script `rounds/027/experimentalist_w.py 8 0 500 600` (575 / 2 / 21 split, same 2 W-and-J != 0 rows). So the two scripts agree on the same walk.
- Not re-run: `enum 7/8`, `reach 7` (2 to 10 min each, 4 jobs); the n = 8 enum counts are not in any saved file (the author says so). The "0 of 42 at n = 6" claim, the pivot of the whole result, rests on `enum 6`, which I did not run either. Weaker than it should be, but the hand table is the actual existence claim and it reproduces.

## True?

Hand table: true as reproduced. Claim (1), existence of a gate-admitted both-die algebra with J != 0 at n = 6, is demonstrated by M1-M4.

Reconciliation of the flags (settled):

1. **136 vs 61 is not a discrepancy.** Both count (algebra, v) rows with J != 0 at out-degree 2 on the n = 8 class 0 walk. E-113's 61 come from a 540 s walk run in parallel (5 114 out-degree 2 rows total); the author's walk expanded 6 855 algebras. On the common deterministic prefix (600 expansions) the two scripts agree exactly. The count scales with walk length, as E-113 itself states ("counts are load-dependent lower bounds"). The author's "different caps, not reconciled" can be replaced by this. Better still: report the count at a fixed `max_exp` (e.g. 6 000) for both scripts.
2. **Out-degree 3 row is a genuine, not-in-E-113 observation, and it is not a conflict.** E-113 says "out-degree >= 3: 3 494 rows, 0 rejects" in a smaller walk. I re-ran the author's walk with a print on out-degree 3 rows (scratch script, not committed): J != 0 at out-degree 3 occurs at expansion 6 820 and 7 424 (v = 4), and 7 822 and 7 831 (v = 7). All four have parallel out-arrows (v -> 8 twice), `mutationIsPossibleAtVertex` True (gate-admitted). So there are at least 4 such (algebra, v) rows in the longer walk, not 1 (the author's count is from a 6 855-expansion walk; mine ran past it). The first two are the author's row. E-113's "0" is a statement about a prefix, so there is no contradiction, but E-113's reading ("W / J != 0 only at out-degree 2") does not extend. Caveat: my mutation check in the scratch script called a non-existent function (`mutation.mutateAtVertex`; it is `procedure.mutateAtVertex`), so I did not confirm that `checkCartan=True` fails on these. The author's caveat (dim count with parallel arrows validated only for n <= 6, E-121, though E-120 line in EXPERIMENTS.md reports a 1 132-quiver check for 4 to 6 vertices) stands for n = 8.
3. The "(3, (1,1))" row has relations with a repeated path ([2,4,8],[2,4,8]); these are coefficient-2 or parallel-arrow terms, not a degenerate relation. Fine, but the author should say that the repeated entry denotes two different parallel arrows.

Gaps:

- Claim (4), "J_i at most 1-dimensional ... always at i with dim e_iAe_v = 2": the n = 8 c0 output shows 179 J_i rows with dim 2, while the claim text says "243 of 243" (c0 + c1 + n = 6, 7?). Check: n=8 c0 179 alone in the file; n = 6 and 7 are all (1) with one i each (234 + 115). 179 + 56 (c1 counted as ~64?) does not obviously make 243. Give the per-file split. Also "dim J_i <= 1" is a statement over capped walks, and is not tied to a mechanism; as stated ("in every case") it is only "in every case found".
- Claim (2) quantifier: "0 of 42 at n = 6" is over the author's enumeration (4 core shapes, one-arrow pendants, scalar 1, monomial kills). The title says "its Coxeter key is not an LNA key at n = 6" without that qualifier. Title should say "in the enumerated families". The hand table does show the key non-LNA for M1-M4, which is solid for those algebras; the generalisation is not.
- The observation that arrow-side cores (B, C) never get an LNA key while (2,2), (2,3) do at n >= 7: the author leaves the mechanism undone. A one-line Coxeter-polynomial reason may exist (the keys M1-M3 are (1,1,-2,-4,-2,1,1): the square contributes a negative middle term the LNA keys at n = 6 cannot reach), but I did not derive it; do not ask for it here.

## New?

- `grep -rn "both-die"` in `research/`: only E-122 (EXPERIMENTS.md line 23, "Both" rows) uses the pattern under the name Both/both-terms-die; the n = 6 existence with J != 0 and the key comparison are not recorded. E-113 (line 120): out-degree >= 3 "0 rejects" (see above). E-123 per author. New as stated.
- Not checked: RETRACTIONS for the key claim; author says not there.

## Evidenced?

- Hand table: yes, specific and reproduced.
- Enumeration counts: counts given, but n = 8 not saved to a file, and the "0 of 21 736" total for B and C at n = 8 is arithmetic on counts that cannot be checked from the repo. Save the n = 8 output.
- BFS reach "contains 0 of 44": stated as capped, correct hedging; "36 877 expanded, 80 680 seen" is the sum of the two classes (23 092 + 13 785 = 36 877; 46 989 + 33 691 = 80 680, consistent).
- Walk J_i tables: rows not totals, labelled so; fine.

## Required for acceptance

1. Restate the 136 vs 61 mismatch as settled: same row definition, different walk length; give counts at a common `max_exp` for both scripts (the 600-expansion prefix gives 2 and 2).
2. Replace "E-113 conflict unverified" with the fact: at least 4 out-degree 3 rows with J != 0 (all with parallel out-arrows, gate-admitted) in the n = 8 c0 walk beyond ~6 800 expansions; E-113's 0 is a prefix statement. Run `procedure.mutateAtVertex(..., checkCartan=True)` on them to say whether they are real rejects.
3. Qualify the title and claim (2) with "in the enumerated families (4 core shapes ...)".
4. Fix or explain "243 of 243" against the per-file counts (n = 8 c0 file shows 179).
5. Save the n = 8 enum output to a file.
