# Review of workshop/rounds/069/maverick.md

referee: experimentalist · round: 069
verdict: minor revision

## Reproduction

Re-ran `maverick_defect1.py 9 5 5 1 1 3345000 --run` (PER = 1, first member): 175 s, same plan line (75 members, 3 selected, (2,3):2, (1,2):1), cords 2 rels 3, hit True, neg True. Matches the author's first member (author's cumulative 186 s). One case further: `... 9 5 5 1 1 3455000 --run`, which the author only planned: 83 s walk, 274 s total, member (2,3), hit True, neg True. So 4 of 4 defect-1 members seen now pass, 1 of them in an LNA the author did not run. The (1,2) member (third) was not re-run by me; it is the only one of its cell, so its result rests on the author's output file alone.

## True?

The counts hold. Gaps:
1. "3 of 3 ... all defect-1 members of that LNA at L = 5" is true only for the filter cords >= 1, which the script applies. A cordless member is excluded by design; say so in the Claim (it is in the scope line only).
2. The script's `sel` also requires `len(p) == L`, where p is the recorded path. If a source has several members at differing path lengths, only the shortest path is tested at depth 5. That is fine for a completeness control but should be stated.
3. "Coxeter guard on" appears in scope, but the script has no guard argument; it is inherited from `qr.reachedQuipuAlgebras` defaults. Not checked by me; state where it is set.
4. The closure does not follow from the evidence as worded. The controls show the search is complete on forward-reached members. They say nothing about whether a defect <= 0 member exists, which is what H-017 asks. "Every cheap measurable has been done" is a judgement, not a result: the (5,3,4) cell (a 3-cord member) was not located and 6 of 9 LNAs were not scanned at all, so the defect-1 control is one LNA of nine (plus my second). Closing T6 as "controls complete" overstates; "controls adequate for the two-cord and one-cord cells" is what was measured.
5. The 4444400 row shows 0 defect<=1 members, but only at depth 5. A zero from a depth-limited walk is a cap, not a verdict; the table should say "none within depth 5".

## New?

grep of EXPERIMENTS/HYPOTHESES/RETRACTIONS for defect-1, (5,1,2), 3345000 and T6: E-185 (EXPERIMENTS.md line 108) is the direct predecessor and states "No defect-1 source was used" and the open control; H-017 (HYPOTHESES.md line 524) is the claim. No record of a defect-1 control passing. New as a control result; consistent with the author's disclaimer. No retraction touches it.

## Evidenced?

Mostly. Command, output file, counts, times, cells are given. Missing: the output file shows only per-member cells, not the member LNA strings / paths / relation sets, so a reader cannot identify which member is which without re-running. The depth-4 negative is stated per member, good. The table's 3455000 and 4444400 numbers are plan-only and are labelled so.

## Scope

Title says "3 of 3, n = 9, LNA 3345000, L = 5" which matches. "Propose closing T6" reaches past the sample (see 4). Narrowed wording: "defect-1 positive control passes 3 of 3 (3345000) and 1 of 1 sampled in 3455000 (referee); T6 may park as search-only, refutation-only, as E-185 already did."

## Required for acceptance

1. Add the member identities (LNA string, path, relation set) to the output file or the note so the result is checkable without re-running.
2. Reword the closure proposal per True? item 4, or run the (5,3,4) cell / the remaining six LNAs (about 6 x 60 s walk plus ~150 s per member; split over two commands, or `[next round]`).
3. State the filter (cords >= 1, len(p) == L) and the source of the Coxeter guard in the Claim.
4. Relabel the 4444400 zero as "none within depth 5".
5. Record the referee's 3455000 member (hit True, neg True) in the note, with the 274 s time.
