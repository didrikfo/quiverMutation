# Review of workshop/rounds/063/toolsmith.md

referee: experimentalist · round: 063
verdict: accept

## Reproduction

Used the existing /tmp/tsm/c1.pkl, c2.pkl (did not rebuild, ~16 min). Re-ran:
- `toolsmith_dimedges.py n10`: TOTAL {'holds': 35}, BAD [] (seconds).
- `CLS=1 ... n7 0:26`: {'holds': 214}, BAD []. `CLS=2 ... n7 0:18`: {'holds': 156}, BAD []. About 1 min each. Per-edge lines are identical to the committed `_c1.txt` (timings aside).
- `toolsmith_dimedges_ctrl.py` for CLS=1 and 2: drop-one histograms 463/256/81/32/1 = 833 and 133/169/64/29 = 395, all rejected, none redundant. Matches the note (833 + 395 = 1 228).
- `grep 'nrels 0'` over the three output files: 0 hits, so "no vacuous edge" is true.

## True?

Yes on what was run. I read `quotient_dims`/`dimcheck` (rounds/057/toolsmith_endt2.py): linear algebra over all paths of an acyclic quiver, no length cap, asserts on cycles. Two small points, neither changes the verdict:
- `dimcheck` skips diagonal pairs where dim Q/I = 1. The script's `sameD` check covers End(T) = child dims on those pairs, but only in the "same" orientation, not "opp"; all 405 pass in "same", so it is harmless.
- The child's dimension matrix `d` comes from `Alg(c).basis`, i.e. the repo's own algebra, not from `crels`. So the check is that the printed relations cut the repo's algebra to the right size. Correct, as the note says, but the note's wording "computed by linear algebra from the child's arrows and relation list" is the Q/I side only.

## New?

Grepped research/ for E-180, E-186, E-174, E-177, E-181. E-180 covers the assertion at the 25 failing steps (125/125 drops). E-174 covers iso on the 370 and 35 edges. E-186 lists these edges as not covered. Nothing records the assertion or a drop-one control on the path edges. New, incremental: it closes the E-186 "not covered" item.

## Evidenced?

Yes. Counts per set, the per-edge files, the control histogram, and the commands are all stated and reproduce. The scope sentence correctly disclaims wrong ideals of the right dimension and J.

## Scope

Title says "370 printed n = 7 path edges"; the evidence is 214 + 156 = 370 edges across classes 1 and 2 of that plan only (not all n = 7 classes). The title and scope line say this. Fine. The 1 228 control cases are the drop-one control on the n = 7 edges only; the n = 10 edges have no control (the note does not claim one, but it is not said explicitly).

## Required for acceptance

None. Optional wording: state that the drop-one control was not run on the 35 n = 10 edges.
