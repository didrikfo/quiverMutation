# At n = 9, one pair: the J = 0 join test finds no join of the certified-inequivalent F-010 pair through depth 6 + 6, a weak specificity datum (equivalent controls need up to 6 + 6 or miss at 4 + 4), and no J != 0 step occurs in the balls, so it cannot test the J = 0 premise

author: experimentalist · round: 057 · kind: negative
thread: T10 · bears on: H-015, F-047, F-010, E-155, E-159, E-161, E-167, E-170
scope: n = 9 only; one certified-inequivalent equal-key pair (3060000 vs 3304000, F-010 / F-047) plus 3 variants; walks of depth <= 6 per side from each LNA and its dual, gate + Coxeter-key guard; Hom test on depth <= 4. n = 10 groups (E-170: 3 profile-separated groups, the Phi^18 group) NOT run.

## Response to referee (theorist, minor revision)

1. Reproduction: done. Commands must be run from the repository root (the script loads round 054's `experimentalist_amerge.py` by a relative path); stated in the Reproduction block.
2. Claim (3): done. "Explained by the key guard" removed (the pair has equal keys, so the guard cannot separate it). Now says only that the J0 and ALL balls coincide through depth 6 and that all sampled edges are J = 0 and Hom-tilting.
3. Claim (1): done. Stated as one pair, and as a non-join at a depth where equivalent controls also fail (2223030 and the second-quipu pair at 4 + 4) or need 6 + 6 (3060000/6000030); weak evidence of specificity. The inequivalent balls at 6 + 6 (4080/2304) are no deeper than the control's own limit.
4. "The test can answer no" softened: the run is circular given the premise (E-166 gives tilting per loop-free step), so it checks the implementation, not the premise. Also noted: the J0 filter removed nothing here, so identical balls only show it is vacuous in this region; the 1551 figure counts DFS-tree edges, not distinct steps.
5. Depth 8 + 8 run: [next round], left as the OVERNIGHT proposal.

## Claim

(1) Specificity: for the F-010 pair (equal Coxeter polynomial, Z-conjugacy profile of F-047 differs, so certified inequivalent by Ladkani Cor 3.15), the relabelling-aware meet in the middle (WL hash + VF2, as E-167) of the J = 0 `tiltingPlus` balls finds 0 joins at depth 2, 3, 4, 5, 6 per side (balls 4080 and 2304 quivers at depth 6). This is one pair, and a weak datum: equivalent controls also fail to join at 4 + 4 (2223030, second quipu) or need 6 + 6 (3060000/6000030), so a non-join at 6 + 6 may reflect depth and not specificity. Given the premise (E-166) the result is also circular: it checks that the implementation does not contradict it.
(2) Sensitivity: certified-equivalent pairs (same quipu, F-047) with different move orbits join: 3060000 ~ 3030000 at total 5 (27 meetings at 5+5), 3060000 ~ 6000030 at total 9 (178 meetings at 6+6, paths in the table). 3060000 vs 2223030 does not join at 4+4 (not tried deeper): a miss is not a verdict.
(3) What the control does NOT show: in every ball above the J = 0 ball and the gate + key-guard ball are identical (same sizes, same joins), and all 1551 gate+key edges up to depth 4 (both sides, dual too) have J = 0 and Hom(T,T[+-1]) = 0 (0 discordances, 0 J != 0). The balls coincide through depth 6, so the J = 0 filter removed nothing here. The control cannot say what a J != 0 step would do. E-149 puts first failures at parent depth 7-8 (n = 7), outside depth 6; I did not reach them. Refuted if: a J = 0 ball of either row, at any depth, meets a ball of the other.
It does not claim the premise (J = 0 steps are tilting) is safe in general, or anything about n = 10.

## Evidence

Joins (total = forward + backward length; row pair, per-side depth):

| pair | status | depth | A ball | B ball | meetings | shortest total |
|---|---|---|---|---|---|---|
| 3060000 / 3304000 | inequiv. (F-010) | 2,3,4,5,6 | 55..4080 | 43..2304 | 0 at every depth | none |
| 3060000 / 3030000 | equiv. | 4, 5 | 549, 1527 | 440, 1103 | 2, 27 | 5 |
| 3060000 / 6000030 | equiv. | 4; 6 | 549; 4080 | 541; 3992 | 0; 178 | 9 |
| 3060000 / 2223030 | equiv. | 4 | 549 | 79 | 0 | none (miss) |
| 3304000 / 2400230 | equiv. (2nd quipu) | 4 | 357 | 146 | 0 | none (miss) |

Sensitivity is partial: of 4 certified-equivalent pairs, 1 joins at 4 + 4, and 3060000/6000030 joins only at 6 + 6 (total 9); 2223030 and the second-quipu pair were not tried deeper. Example winning path 3060000 ~ 6000030: forward [1,1,3,1,2,2], backward from the target [-8,-9,-9], meeting equal up to relabelling (inverse of the back path: [3,3,2]); 178 meetings, 0 J != 0 steps on the best one. Hom table (`experimentalist_powerhom_d4.txt`): rows 3060000 (right 457, dual 475 edges) and 3304000 (315, 304): all J0 and all tilt. Counts are edges of the DFS tree with memo, not distinct steps.
The asymmetry is the power statement: a join between inequivalent classes would need an edge that is not an equivalence; the Hom test shows none among 1551, and E-166 proves generation per loop-free step, so impossibility of the false join follows from the premise at depth 6, and the run only confirms the implementation does not contradict it.

## Reproduction

```
# run all from the repository root (relative path to round 054's amerge script)
timeout 10m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py 9 3060000 3304000 6        # 310 s (J0 + ALL, 155 s each)
timeout 10m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py 9 3060000 3030000 5 J0     # 51 s
timeout 10m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py 9 3060000 6000030 6 J0     # 175 s
.venv/bin/python workshop/rounds/057/experimentalist_powerhom.py 9 3060000 3304000 4                      # 37 s
```
Outputs: `experimentalist_powerjoin_{d4,d5,d6,ctrl_d5,ctrl2_d6}.txt`, `experimentalist_powerhom_d4.txt` (same dir).

## Prior record

E-170/E-168: certified inequivalent key-equal pairs exist only as the F-010 pair (n = 9) and 3 profile-separated groups at n = 10; the 13 n = 10 equal-profile groups are uncertified; nothing in `research/` runs the join test on such a pair. E-155/E-161 controls (toolsmith_tiltpath ctrlN:L) are sensitivity controls only (random walks of known length), E-167 says its test "is not shown able to fail" (no J != 0 step). This round gives the first specificity datum, and shows the same weakness persists: no J != 0 step in the balls. Nothing in RETRACTIONS touched (grep for "specificity", "power" found none).

## Code changed

None to `quivermutation/`. New: `experimentalist_powerjoin.py` (imports `experimentalist_amerge.py` of round 054 for the WL/VF2 join), `experimentalist_powerhom.py` (execs `rounds/050/skeptic_tilt.py`). No tests run (no library file touched).

## Next

- Why a better control is hard: an inequivalent pair is only reachable through a J != 0 step if the walk goes >= 7-8 deep; at n = 9 that is ball sizes ~10^4-10^5 per side. Proposal for OVERNIGHT.md: F-010 pair, depth 8 + 8, J0 and ALL modes, report first J != 0 edge and Hom verdict on it (est. several hours; depth 6 took 155 s, growth ~2.7x per level gives ~19 min at 7, ~1.5 h at 8, per mode per row).
- The meaningful power test of the J = 0 premise is Hom(T,T[-1]) on the J != 0 edges themselves (E-166: nonzero on all 16 tried); that needs a pair whose ball contains J != 0 steps, e.g. n = 7 c1 parents at depth 7-8, with the other side an inequivalent class: toolsmith.
- maverick: supply the n = 10 rows for the 3 profile-separated groups and the Phi^18 group so the same script runs there.
