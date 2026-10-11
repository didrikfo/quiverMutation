# At n = 10, one certified-inequivalent pair of the Phi^18 group (90000000 vs 50505000): the J = 0 join test finds no join through depth 6 + 6, while a 7-step equivalent pair (34504030 ~ 50505000) joins at 4+4 in the same mode

author: maverick · round: 058 · kind: negative
thread: T10 · bears on: H-015, F-047, F-010, E-172, E-175, E-169
scope: n = 10, ONE inequivalent pair (classes {1,2} vs 3 of the one live Phi^18 key group; 50505000 is derived equivalent to 34504030, F-037/E-032); J0 mode to depth 6 per side, ALL mode (gate + key only) to depth 5; 4 equivalent control pairs at depth <= 4, one of them a 7-step path; no Hom test, no other pair, no depth >= 7.

## Response to referee

1. Control 34504030 ~ 50505000 (4+4, J0): done, 31 s, reach 275/588, 28 meetings, shortest total 7, 0 J != 0 steps on both paths (`maverick_join_equiv_long.txt`, row added to the table). Statement "a control ... is not available" deleted. Consequence: the J0 search at 4+4 finds a 7-step equivalence from 50505000's side, so the empty 90000000/50505000 join at 6+6 (B ball 5438 > 588) is stronger evidence than the first draft said. Caveat: same B ball as the test, but the test's A side (90000000) has no equivalent control; a control for that side is still missing. The toolsmith item is restated for the 90000000 side only.
2. Claim (1) corrected: 34504030 ~ 50505000 (F-037, E-032, H-013), so there are 3 derived-equivalence candidates {0}, {1,2}, {3} and 2 distinct certified inequivalences ({0} vs {1,2}, {1,2} vs {3}); 0 vs 3 share a profile (not certified). "4 certified pairs" and "singleton" wording removed (only F-047-move singletons).
3. Prior record sentence corrected (F-037/E-032/H-013 hits for 50505000; 90000000 not mentioned in research/).
4. J0 ball vs gate+key ball: ALL mode run at depths 3, 4, 5 (`maverick_join_ineq_ALL_d3to5.txt`): reach 373/185, 1325/588, 4273/1823, identical to J0 and 0 meetings. So the J0 ball equals the gate+key ball through depth 5: no J != 0 step enters these balls, and J0 and ALL are the same test here. Depth 6 ALL not run (about 7 min; same trend expected, not verified). The "30000002 misses" is no longer the best calibration.

## Claim

(1) The E-172 live group (Coxeter polynomial (T+1)^2(T^2-T+1)(T^6-T^3+1), Phi^18 = I) has 4 classes (orbit + mirror of F-047 moves): class 0 = 320 rows (one orbit, includes 00000030), and three singleton classes 34504030, 50505000, 90000000. Computed with the F-047 Smith profile: the profile sets of 0 vs 1, 0 vs 2, 1 vs 3, 2 vs 3 are disjoint (certified inequivalent by the profile criterion used for F-010, Ladkani Cor 3.15); 0 vs 3 and 1 vs 2 share one profile (not certified). Since 34504030 ~ 50505000 (F-037/E-032), this is 2 distinct certified inequivalences among the candidates {0}, {1,2}, {3}; the group is not "4 mutually certified classes".
(2) For the pair 90000000 (class 3) vs 50505000 (class 2), the relabelling-aware J = 0 meet (as E-175) has 0 joins at depth 3, 4, 5, 6 per side (reach sets 373/185, 1325/588, 4273/1823, 12855/5438 quivers).
(3) Sensitivity at n = 10 is only partial: inside class 0 (one orbit), 00000030 ~ 00000230 joins at total 4 (66 meetings at 4+4), 00000030 ~ 20000030 joins at total 1, but 00000030 vs 30000002 (same orbit) does NOT join at 4+4. A miss at depth 4 is not a verdict, so the depth-6 non-join of (2) is weak evidence of specificity, no more than E-175's.
It does not claim: the J = 0 premise holds; the pair is inequivalent beyond the profile criterion; any statement at depth >= 7, where E-151 puts the first J != 0 failures (at n = 7). I did not count J != 0 steps (J0 mode prunes them), so I cannot say whether any lies in these balls. Refuted if: any J = 0 ball of 90000000 meets one of 50505000 at any depth.
Power size: the balls grow about 3.0x per level (4273 -> 12855 and 140 s -> 424 s on the A+B pair). A plausible join test, for a pair whose controls need 6+6 at n = 9 and need a J != 0 step, is depth 8+8: about 9x the depth-6 cost, ~1 h per mode per pair (single core), more if the dual walk memory grows. Not run.

## Evidence

| pair | status | depth | reach A / B | meetings | shortest total |
|---|---|---|---|---|---|
| 90000000 / 50505000 | profile-disjoint | 3 | 373 / 185 | 0 | none |
| same | | 4 | 1325 / 588 | 0 | none |
| same | | 5 | 4273 / 1823 | 0 | none (140 s) |
| same | | 6 | 12855 / 5438 | 0 | none (424 s) |
| 34504030 / 50505000 | derived equivalent (F-037), 7 steps | 4 | 275 / 588 | 28 | 7 (31 s) |
| 90000000 / 50505000 | ALL mode | 3, 4, 5 | same as J0 | 0 | none |
| 00000030 / 00000230 | same orbit | 3, 4 | 351/247, 1194/750 | 13, 66 | 4 |
| 00000030 / 20000030 | same orbit | 3, 4 | 351/271, 1194/871 | 98, 487 | 1 |
| 00000030 / 30000002 | same orbit | 3, 4 | 351/271, 1194/871 | 0, 0 | none (miss) |

Control paths are one-to-three steps: the 3-in-class controls test only that the script runs at n = 10 and joins close pairs; the pair that needs a long path (30000002) misses at 4+4. 
Reading: with no J != 0 data, the empty join at 6+6 is what the premise predicts and is also what a weak search predicts. Correct control for the next level: a pair certified equivalent that needs a long path (n = 9: 3060000 ~ 6000030 needed 6+6; n = 10 analogue: 34504030 ~ 50505000, total 7, joins at 4+4).

## Reproduction

```
.venv/bin/python workshop/rounds/058/maverick_group.py                         # 17 s: classes, profile disjointness
.venv/bin/python workshop/rounds/058/maverick_ctrlpairs.py                     # orbit of 00000030 (320 rows)
timeout 10m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py 10 90000000 50505000 6 J0   # 424 s (also depth 3, 4, 5: 11 s, 39 s, 140 s)
timeout 10m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py 10 00000030 00000230 4 J0   # 32 s (controls: also 30000002, 20000030)
```
Outputs: `maverick_join_equiv_long.txt`, `maverick_join_ineq_ALL_d3to5.txt`, `workshop/rounds/058/maverick_join_ineq.txt`, `maverick_join_ineq_d6.txt`, `maverick_join_ctrl.txt`.

## Prior record

E-175 (n = 9, F-010 pair, 6+6, J0 and ALL), E-169 (n = 10 join machinery), E-172 (the group: char poly, Phi^18, 4 classes; it did not list rows or which pairs the profile certifies). Grep of research/: 50505000 appears with 34504030 (F-037, H-013, E-032, bruestle/2310.08346 notes: derived equivalent by seven mutations); 90000000 nowhere. The n = 10 certified-inequivalent pairs of E-172's other three profile-separated groups are untouched.

## Code changed

None to `quivermutation/`. New: `maverick_group.py`, `maverick_ctrlpairs.py`. The join script of round 057 runs unchanged at n = 10. No tests run.

## Next

- experimentalist: depth 7+7 and 8+8 for this pair in J0 and ALL modes as an OVERNIGHT job (~1 h and ~9 h projected single-core; a shard by side is possible), reporting the first J != 0 edge.
- toolsmith: a long-path equivalent control for the 90000000 side only (50505000's side is covered by 34504030).
- The 90000000 class (one relation of length 9 on A_10) and 50505000 being singleton classes suggests a direct proof route (their tilting complexes are very restricted); theorist.
