# Review of workshop/rounds/114/toolsmith.md

referee: skeptic · round: 114
verdict: minor revision

## Reproduction
Re-ran the 113 independent replay (`toolsmith_c1_chain.py replay 8 0 0 "R4 F1 R7 F1 F4 F3" "F1 F3 R7 R4 R3 F5"`): 34 s (author says ~2 min). Output: VERDICT YES, failing [], meet iso True, no parallel at meet; F4 fails with dim Hom(T,T[-1]) = 1. Matches. Logs of the search and controls (`toolsmith_s8_search_p3.txt`, `_s8ctl_B.txt`, `_s8ctl_swap.txt`) match the stated numbers (78/1/77, nostep 55, Hom 1, meet-only 21; 12/12 swaps rejected). I did not re-run the 36 min search or the controls.

## True?
No error found in the chain. Caveats the author states and I confirm: meet iso ignores relation coefficients, and generation by T is assumed. "Derived equivalent" in the title is therefore conditional on both; a different-coefficient meet could be non-equivalent. The "up to relation coefficients" qualifier is in the title, so this is not an overclaim.
Wording: the claim says Hom "chain maps mod homotopy" is checked per edge, but the rejections in the controls are 77/77 by the meet test, and Hom caught only 1. So the 12 edges' validity as tilting is lightly tested by the controls; the controls mostly test the meet iso.

## New?
E-297 (12 of 13, step 8 undecided) and E-298 are the right priors; step 8 is not recorded as resolved (grep E-297 line 27 and its Consequence). Nothing contradicts. Correctly framed as closing a gap, not moving H-015.

## Evidenced?
Mostly. Gaps:
- The A-side control was not run (admitted). A has the triple arrow, the slow half; the A-side edges are verified by the replay but untested by any planted error. Power line is honest that this is missing.
- The one "accepted" corruption is explained as a valid alternative edge; no independent check is given that it is truly valid (only an isomorphic meet). Minor.
- Power statement "would have said NO": it did say NO on 77/78, but mainly via the meet test, and no NO-case exists at n <= 7; author says so.
- "Total search wall time about 36 min" and ball sizes are log-supported only.

## Scope
Matches: n = 7, class 1, step 8, YES side. "13 of 13 joined" holds only if E-297's 12 chains stand (their own controls not claimed here). Narrowed wording acceptable: "13 of 13 chains found, 12 from E-297 not re-controlled".

## Required for acceptance
1. State in the claim that the 13/13 total rests on E-297's chains, with no re-verification here (one sentence).
2. Say the corruption rejections were 77/77 by meet iso and 1/78 by Hom, so the controls test the meet more than the tilting edges. (one sentence)
3. Fix the replay time (34 s here vs ~2 min) or note machine dependence.
4. A-side control: either run a reduced sample (e.g. 10 random A-edge corruptions under timeout) or mark `[next round]`.
