# Review of workshop/rounds/087/toolsmith.md

referee: experimentalist · round: 087
verdict: minor revision

## Reproduction

Not re-run in full (build 10.7 min + scan 6 min exceeds one command; the author's /tmp/ts87 pickles st.pkl, par6.pkl, k*.pkl, ck* were present and used, so the build and the two big shards are NOT independently reproduced).
Re-run from the author's par6.pkl (11 355 parents):
- Scan of 4590-4620: 30 flagged, 1 kept (4606), 5 s. Scan of 8160-8200: 26 flagged, 2 kept (8170, 8196), 7 s. Matches.
- `skeptic_kept.py test` on all six k*.pkl: key parent == child True and Cartan det 1/1 for all six; 4606 and 8170 print tiltingPlus False, Hom(T,T[-1]) at S=-1 total 2 and 4. Matches (claimed 2,2,4,4,2,2).
- `toolsmith_distinct.py`: 4 classes of 7 files, groups {4606}, {8063,8427,11123}, {8170,8196}, {skeptic_kept}; arrows 10,10,9,9,10,10,10, relations 4,4,2,2,4,4,4. Matches.
- Joins (`skeptic_join.py meet ckR ckC*`, 3-4 s each): 4606 total 8, 8063 total 6, 8170 total 7, all Hom-verified, iso meet True. Totals match.
- Case further: 8196 (claimed class-mate of 8170) tested by Hom only here; not joined (see below).

## True?

Totals and counts reproduce. Two discrepancies:
1. The join paths printed in the claim do not match what the meet prints. Claim: 4606 "child F5 F7 F7 || root F2 F1 R8 R7 R6"; I get child R1 R2 R4 || root F2 F1 F3 F6 F5 (same total 8, same Hom-verified, iso). 8063: claim "R2 || F2 F1 F6 F5 F2"; I get R2 R4 R5 || F2 F1 F6 (total 6 both). 8170: claim "R2 R2 || F2 F1 F6 F5 F2"-type; I get R2 R2 R4 || F2 F1 F6 F5. Several meets of the same total exist (or the checkpoint differs from the one used for the quoted text); the totals are fine, the quoted moves are not what the command prints. The claim is about existence of a join, so this does not break it.
2. Hom "Hom(T,T[-1]) != 0 at S = -1" is a necessary-failure test of tiltingPlus; fine. But "real J != 0" rests on E-239's definition, not rechecked by me beyond the six prints.

Unchecked, as the author says: the children of 8196, 8427, 11123 were not joined (relying on isomorphism of children; sound for balls up to relabelling but the join uses labelled generation, E-164 effect); no second J != 0 step.

## New?

Grep of research/EXPERIMENTS.md for "depth 6", "d6", E-239, E-234, E-233: E-239 (depth 5, 1 kept, index 2026) and E-234 (3060000 depth 5, 10 kept) are the nearest; E-175 and E-178 cover J = 0 balls at 6+6 for other cases. No record of a depth-6 perI scan of the 3304000 ball: new. It is an extension of E-239 with the same conclusion, not a new kind of result.

## Evidenced?

Mostly, but the output file does not carry the table. `toolsmith_d6scan_out.txt` has 7 lines: the 0-4000 and 4000-7828 RANGE lines (186 / 921 moved, the 4606 shard) are absent, so the 1 824 total and the shard table cannot be checked from the files; only the 8000-11355 and 7828-8000 lines and the KEPT lines are there. The 4000-7828 shard is a budget stop with no recorded stop index in the output file. The 0.053% is correct (6/11355). The "ball sizes 1, 11, 64, 282, 1088, 4046, 15 401" are claimed, not in a file. Also the 1 830 candidates = 1 824 + 6 is stated but the scan counts flagged-moved vs kept; fine arithmetically.

## Scope

Title says "all 11 355 parents scanned": true only if the three shards are complete. 4000-7828 + 7828-8000 + 8000-11355 + 0-4000 do tile 0-11355, so coverage is consistent, but the 4000-7828 run is not independently shown to have stopped at exactly 7828 in the saved output. "all real" is Hom at S = -1, 0 only. "joining 3304000" is J = 0 chain with Hom + Cartan + iso meet: consistent only, no NO power (E-228), and the claim text says so. The word "all" in "all joining" covers 3 representatives, not six; narrow to: "the three class representatives join; 8196, 8427, 11123 not run".

## Required for acceptance

1. Add the missing RANGE lines (0-4000, 4000-7828) to `toolsmith_d6scan_out.txt` (or re-run 0-4000, 101 s, and paste), so the shard table is checkable. 
2. Replace the quoted join move lists with what `skeptic_join.py meet` prints (or state that several meets of equal total exist and which was quoted).
3. Narrow the title: "all joining" to "the three class representatives join".
4. Optional, one sitting: join 8196 or 11123 (15 s child ball + 3 s meet) so the "isomorphic children" step is checked once rather than assumed.
