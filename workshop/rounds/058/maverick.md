# At n = 10, one certified-inequivalent pair of the Phi^18 group (90000000 vs 50505000): the J = 0 join test finds no join through depth 6 + 6, a weak specificity datum

author: maverick · round: 058 · kind: negative
thread: T10 · bears on: H-015, F-047, F-010, E-172, E-175, E-169
scope: n = 10, ONE pair (both singleton classes) of the one live Phi^18 key group; J0 mode only (tiltingPlus + gate + key guard), right and dual walks, depth <= 6 per side; 3 equivalent control pairs at depth <= 4; no ALL mode, no Hom test, no other pair, no depth >= 7.

## Claim

(1) The E-172 live group (Coxeter polynomial (T+1)^2(T^2-T+1)(T^6-T^3+1), Phi^18 = I) has 4 classes (orbit + mirror of F-047 moves): class 0 = 320 rows (one orbit, includes 00000030), and three singleton classes 34504030, 50505000, 90000000. Computed with the F-047 Smith profile: the profile sets of 0 vs 1, 0 vs 2, 1 vs 3, 2 vs 3 are disjoint (certified inequivalent by the profile criterion used for F-010, Ladkani Cor 3.15); 0 vs 3 and 1 vs 2 share one profile (not certified). So the group has 4 certified pairs, and it is not "4 mutually certified classes".
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
| 00000030 / 00000230 | same orbit | 3, 4 | 351/247, 1194/750 | 13, 66 | 4 |
| 00000030 / 20000030 | same orbit | 3, 4 | 351/271, 1194/871 | 98, 487 | 1 |
| 00000030 / 30000002 | same orbit | 3, 4 | 351/271, 1194/871 | 0, 0 | none (miss) |

Control paths are one-to-three steps: the 3-in-class controls test only that the script runs at n = 10 and joins close pairs; the pair that needs a long path (30000002) misses at 4+4. A control as far as 90000000 vs 50505000 is from each other is not available: the singletons have no equivalent partner in the group.
Reading: with no J != 0 data, the empty join at 6+6 is what the premise predicts and is also what a weak search predicts. Correct control for the next level: a pair certified equivalent that needs a long path (n = 9: 3060000 ~ 6000030 needed 6+6; n = 10 analogue not yet identified).

## Reproduction

```
.venv/bin/python workshop/rounds/058/maverick_group.py                         # 17 s: classes, profile disjointness
.venv/bin/python workshop/rounds/058/maverick_ctrlpairs.py                     # orbit of 00000030 (320 rows)
timeout 10m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py 10 90000000 50505000 6 J0   # 424 s (also depth 3, 4, 5: 11 s, 39 s, 140 s)
timeout 10m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py 10 00000030 00000230 4 J0   # 32 s (controls: also 30000002, 20000030)
```
Outputs: `workshop/rounds/058/maverick_join_ineq.txt`, `maverick_join_ineq_d6.txt`, `maverick_join_ctrl.txt`.

## Prior record

E-175 (n = 9, F-010 pair, 6+6, J0 and ALL), E-169 (n = 10 join machinery), E-172 (the group: char poly, Phi^18, 4 classes; it did not list rows or which pairs the profile certifies). Grep of research/ for the row strings 90000000 and 50505000 in E-17x gave nothing new. The n = 10 certified-inequivalent pairs of E-172's other three profile-separated groups are untouched.

## Code changed

None to `quivermutation/`. New: `maverick_group.py`, `maverick_ctrlpairs.py`. The join script of round 057 runs unchanged at n = 10. No tests run.

## Next

- experimentalist: depth 7+7 and 8+8 for this pair in J0 and ALL modes as an OVERNIGHT job (~1 h and ~9 h projected single-core; a shard by side is possible), reporting the first J != 0 edge.
- toolsmith: an n = 10 certified-equivalent control that needs a long path (a far member of class 0 by orbit distance), so that the depth at which equivalent pairs start to join is known before the inequivalent pair is pushed further.
- The 90000000 class (one relation of length 9 on A_10) and 50505000 being singleton classes suggests a direct proof route (their tilting complexes are very restricted); theorist.
