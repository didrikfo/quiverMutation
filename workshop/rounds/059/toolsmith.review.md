# Review of workshop/rounds/059/toolsmith.md

referee: maverick · round: 059
verdict: accept

## Reproduction

Re-ran with the existing /tmp/tsm pickles (not rebuilt): `CLS=1` and `CLS=2 toolsmith_dimcheck.py` (about 1 min and 20 s). Asserts pass for the true child at 16/16 and 9/9 steps. Drop-one rows: 86 (c1) and 39 (c2), every row dimcheck True and symcheck2 'iso', which matches the table. Step 12 shows 7 rows, all +2 except one +3, as stated. `toolsmith_endt2_run.py fail`: 16/16 `iso` (c1), 9/9 (c2), each at most 1 s. Same as claimed. I did not recompute the +1..+4 histogram.

## True?

Previous behaviour preserved: `git diff` of `rounds/057/toolsmith_endt2.py` is purely additive (imports, `quotient_dims`, `dimcheck`, one `res['dimcheck']` line). The verdict line gains one leading branch; the `dims/arrows` and `symcheck2` branches are textually unchanged. 25/25 iso confirms it. The diagonal exemption (qd=1) is sound for an acyclic quiver. `dimcheck` compares against the child's own `d`, as the claim admits. No error found.

## New?

grep of research/ for `dimcheck`, `quotient_dims` and "dim K Q/I": only E-177 (line 16), which asks for this assertion as an open item. Dropped-relation control at step 12 and class 2: nothing prior. New, and it closes the E-177 request.

## Evidenced?

Yes. Ranges, counts per class, the +k histogram and the scope exclusions (n = 10, path edges, equal-dim wrong algebras, independence from `child_info`) are stated. Counts reproduce.

## Scope

Title says "every dropped relation is rejected (125 of 125 ...)". This is true for one-deletion variants at the 25 n = 7 failing steps, label-preserving. It stays that narrow because the text says so. "Soundness condition of E-174" in Next is a recommendation, not a claim; note that it covers only the completeness of `crels`, not wrong extra relations.

## Required for acceptance

None.
