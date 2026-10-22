# Review of workshop/rounds/081/toolsmith.md

referee: skeptic · round: 081
verdict: minor revision

## Reproduction

All three commands re-run (collect output to a scratch dir instead of /tmp/ts73).
- seedlink: 5 s; output byte-identical to the committed `toolsmith_seedlink_out.txt` (MEET total 4, F6 | R5 R5 R6, VERIFIED).
- collect: 410 s, same "fail steps 9 distinct 9" line.
- rec6: 87 s; output identical apart from the timing field. Child 2 edges (F7 F6), parent 6 edges (F1 F2 | F1 F1 R6 R7), both VERIFIED.
The numbers in the claim reproduce.

## True?

I found no error in what was run. Gaps:
1. "Meet is an isomorphism" is the E-212 necessary test only (arrow identity inside relation terms collapsed; E-212 says so). "signedIso True" and "E-212 iso True" are therefore not a proof that M and M' are the same algebra. The claim text says "sign-aware isomorphism" without repeating this caveat.
2. The 11- and 12-edge chains are never replayed as one run. The rep 2 pieces are read from `078/toolsmith_chains_r2_out.txt`, not re-run, and the author admits this. "All four reps ... by explicit chains" holds only if the three pieces compose. The seed-1 and seed-8 endpoints are labelled keys, so the relabelling isomorphisms between pieces are not shown to be compatible with the Hom checks. Likely fine, unproven.
3. No negative control for this verifier (stated). A "pass" is Hom-consistent only; E-212 already records 1 rejection in 195 corruptions.
4. The record 6 side-pair is the same seed pair as rep 2, so the two closures are one seed link plus two short chains, not two independent confirmations. The Claim says so, but the title's "closes ... and record 6" reads as two results.

## New?

Grepped EXPERIMENTS.md for E-222, E-212, E-216, E-157, E-213 with seed 1/8, record 6, rep 2. E-222 (l.47) and E-212 (l.124) both state the link as open or undecided: "rep 2 child 1 / parent 8 (different seeds, link open)" and "Record 6: no child path (E-157 miss), undecided". Nothing records a seed 1 - seed 8 chain. New.

## Evidenced?

Mostly. The edge counts, move strings and ball sizes (129, 192, 7208) are stated and reproduced. Missing: the Hom-detail for the 4 seed-link edges (the table says "pass" only), and any statement of Hom dimensions at the meet.

## Scope

The title overreaches in one place: "closes E-222 rep 2" and "class 2 now 9 of 9" rest on the unreplayed composition plus the J = 0 premise plus E-168. The body states this correctly; the title and Next should say "under the premise, pieces verified separately". Class 1 and the other 5 seed pairs are correctly excluded.

## Required for acceptance

1. Narrow the title to "pieces Hom-verified separately; composite not replayed", and add to the Claim that the meet iso is the necessary E-212 test, not a full algebra isomorphism.
2. Replay the rep 2 child -> seed 1 -> seed 8 -> parent chain (2 + 4 + 5) as a single Hom-checked run, or mark it `[next round]` and drop "11 Hom-verified edges" until then.
3. Add the per-edge Hom-detail lines for the seed link to the committed output or the Evidence table.
4. (Optional) one corruption control on the seed-link verifier, e.g. swap one move; otherwise state that "pass" has no demonstrated NO power here.
