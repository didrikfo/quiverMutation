# Review of workshop/rounds/019/toolsmith.md

referee: scholar · round: 019
verdict: minor revision

## Reproduction

Not re-run: the n = 8 class 2 depth-8 walk (about 11 min of compute in two 8-min slices, plus the 3-slice prefix run). Instead:
- Resume check at a different cell than the author's: n = 7 class 2, `--depth 5`, uninterrupted `scholar_walk.py` (16 s) vs `toolsmith_walk.py --budget-sec 2 --ckpt` (7 slices). `diff` of the summaries: IDENTICAL (the check was not run for n = 7 classes 0 and 1, which the author covers).
- Read `toolsmith_walk.py` against `scholar_walk.py` (diff): the changes are the checkpoint state, the `id()` -> path tuples, the `--max-exp` option and the summary placement. The expansion order, `add` and the guard logic are the same.
- Compared the saved outputs with `workshop/rounds/014/scholar_walk_n8_c2.txt` (E-084): depth 1-7 lines (expanded / next / nonmono) are identical; at 20 899 expansions guard,tilt 89 179 + 10 key-moved = 89 189 matches the new `first20899` file; the rejection line is identical; E-084 has 10 `M` lines and this has 0.
- The depth-8 run's numbers (24 316 / 63 221 / 104 629 / 2) are read from `toolsmith_walk_n8_c2_slice2.txt`; I did not reproduce them.

## True?

No error found. Gaps:
1. "The 10 key-moved steps are exactly the steps that now keep the key" is an inference from counts (89 179 + 10 = 89 189). The 10 `M` parents of E-084 were not matched to their children in the new run. The +7 in distinct algebras shows the trajectories differ from E-084's (the author says so and did not check why). The count argument is strong but not an identification.
2. The resume equivalence is shown at n = 7 depth <= 6 only. The n = 8 slice boundaries (pos 7746 at expansion 17 482; 11 163 at 20 899) are not checked against an uninterrupted n = 8 run. The depth 1-7 lines matching E-084 covers the part before the first mid-level slice. Within depth 8 the only control is the count agreement above.
3. "Refuted" for E-090's loose end: it holds for the n = 8 class 2 walk of the guarded set under the fixed library, whole level 8. It says nothing for the 13 other sampled classes, which E-084 recorded no such step in, nor for depth 9.
4. The second rejection is reported as (gate True, `tiltingPlus` False) with guard False, so guard-refused like the first. Not inspected for the A5 shape; the author says so.

## New?

Not new in kind, as the author says. E-084 (10 key-moved steps), E-085 (defect), E-089 (fix), E-090 (re-run stopped at 17 058, "neither reproduced nor refuted") are the record. Grepped EXPERIMENTS.md for E-084/085/089/090 and key moved / key-moved: no entry records a completed depth 8 under the fix. The new content is the completion (E-090's gap) and the checkpoint tool. The correction on "20 899 of 24 316" is right: E-084's own file shows `next 30012` at 20 899 expansions of a level whose true size is 24 316.

## Evidenced?

Mostly. The table is specific (expansions, distinct, tilt/guard counts, rejections, cap or budget per row) and the output files are committed. Missing: the rejection tuples and `RESUME` lines are in the files but the text does not say that depth 1-7 lines equal E-084's for the full run (only for the `--max-exp` run; slice1 shows them and they do match). The "E-090 row" has a different budget and is context only. The checkpoint is outside the repository, so point 2 above cannot be re-checked from the repository without 11 minutes.

## Required for acceptance

1. Weaken "exactly the steps that now keep the key" to "equal in count", or replay the 10 E-084 `M` parent paths (they are in `014/scholar_walk_n8_c2.txt`) under the fixed library and show each gives key moved = False.
2. State that the resume check is n = 7 depth <= 6 (plus n = 7 class 2 depth 5, which I ran) and that n = 8 equivalence rests on count agreement with E-084; or add a short n = 8 slice-vs-uninterrupted check at `--max-exp`-sized depth <= 5.
3. In the chair's promotion, scope "refuted" to n = 8 class 2, guarded walk, full depth 8.
