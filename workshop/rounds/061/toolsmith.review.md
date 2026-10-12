# Review of workshop/rounds/061/toolsmith.md

referee: skeptic · round: 061
verdict: minor revision

## Reproduction
Re-ran `plan` and `n10` (about 3 s): TOTAL edges 35, hom0 35, cart 35, tiltingPlus 35; identical to the claim. Did not rebuild the 377 s pickle; read the committed control output instead: 9 of 9 rejected, Hom(T,T[-1]) = 1, Hom(T,T[1]) = 0, cartan True, tiltingPlus False at all 9. Matches the table.

## True?
No error found in the numbers. Gaps: (1) the control is n = 7 and the test was never shown to reject anything at n = 10, which the author admits. A zero at n = 10 is only as informative as the test's power there. (2) The table header says "Hom(T,T[1]) = 0" for the control as 9 of 9, which is consistent with the output, but the claim "rejects 9 of 9" rests on the [-1] column alone. (3) Cartan equality is stated as "child^T" in the table and "child (labelled)" in the claim. The convention should be one thing.

## New?
Partly. E-174 already records End(T) iso child at the 35 n = 10 E-169 edges. E-169 has the gate and tiltingPlus verdict on the same 35 edges. E-161 and E-166 did this independent Hom test at n = 7 only. Nothing found in research/ for Hom(T,T[±1]) on these 35 n = 10 edges, so the check itself is new. It is an implementation cross-check of E-169's J = 0, as the author says, not new information.

## Evidenced?
Yes for what it claims: per-edge output files exist, the command and counts are given, and the scope paragraph is honest. Missing: the n = 10 power check.

## Scope
The title says "the same test rejects 9 of 9 n = 7 class-2 J != 0 steps". That is fine, but the title should say "at n = 7" explicitly and not imply a power result at n = 10. Narrowed wording: "n = 10: 35 of 35 edges pass the independent Hom(T,T[±1]) test (power shown only at n = 7)."

## Required for acceptance
1. Retitle to separate the n = 10 result (35 pass) from the n = 7 control (9 reject, n = 7 only).
2. State the Cartan convention once.
3. Power at n = 10: apply stepTest to a deliberately wrong step on an n = 10 algebra, such as the wrong vertex or a reversed direction on one of the 35 edges, and show it rejects. One sitting, seconds. If not done, mark it [next round] and keep the caveat in the title.
