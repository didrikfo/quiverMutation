# At n = 9 (F-010 key group), all 19 recorded key-kept perI J != 0 steps (E-233, E-234, E-239, E-244, E-254) keep q_S, and every child was already joined to its parent's J = 0 component (E-234, E-239, E-244, E-254, E-255, E-268), so q_S could have said NO only if H-023 failed: no recorded step has an outcome H-023 does not predict; E-306's balls hold 2 of them (1 + 1 + 0, as E-306 says)

author: skeptic · round: 121 · kind: negative (mostly rediscovery of E-306's count and of E-254/E-268's component membership, plus the q_S test on all 19 recorded steps)
thread: T10 · bears on: H-015, H-032, H-023
scope: n = 9, F-010 key group only; E-212 `moves` balls: 3060000 depth 4, 3304000 depth 5, 2400230 depth 3 (parents of depth <= ball depth), plus the 17 further recorded steps replayed from their paths (E-234: 10 at 3060000 depth-5 parents; E-244: 6 at 3304000 depth-6 parents; E-254: 1 at a 2400230 depth-8 parent); only perI-flagged one-vertex F/R J != 0 candidates (multi-vertex J != 0 steps are not seen); 3060000 depth-5 ball and 2400230 depth-8 ball not rebuilt (their steps replayed from paths), 3060000 depth 6 not run.
power: on the real steps none. The test fires only if a child leaves its parent's J = 0 component and H-023 fails there; no recorded child leaves it. The planted cross comparison (child vs the other half's root invariant, NO 19 of 19) is E-306's single half-separation bit, INV(3060000) != INV(3304000), repeated, not 19 independent controls; it does not show that the parent/child step test can fire.

## Response to referee (experimentalist, minor revision)

1. **E-254 missing (18 -> 19).** Done. Added E-254 (2400230, path `F1 F1 F2 F3 R9 R7 R9 F3`, step F3) to STEPS and re-ran `replay` (19 s): path rebuilds with `moves` moves, tag `kept`, parent and child in half B, q_S equal, planted NO True; `TALLY {('kept','YES (q_S equal)',True): 19}`. 18 -> 19 in title, scope, P3 (registered on 18; the 19th is a post-registration addition, not a hit), table and Claim. Thank you: I grepped E-233/234/239/244 and not E-254, the same lesson as 8427/E-255.
2. **"R at v twice" shape.** Done, wording corrected and the shape now checked in code. `replay` prints a shape key per step: does the step repeat a move (same kind, same vertex) of its recorded path? Result: 17 of 19 repeat the path's **last** move exactly (16 R after the same R, E-254 F3 after F3); 2 do not occur in the path at all (E-234 4722: `... F1 F8`, step R1; E-244 11123: `... F2 F7`, step R2). So "all of the R-at-v-twice diamond shape" was false. The key is a property of the recorded (shortest, BFS-order) path, not of the step, so I do not call it a shape of the steps; removed from (c) and from the theorist question in Next.
3. **"Predetermined" / "implied before it was computed".** Done. Reworded to "predicted by H-023 (OPEN) and E-306's record". Agreed: E-306's 14 815 constant-q_S steps are ball steps, not the join-chain steps through the child balls.
4. **Power line.** Done, rewritten as the referee suggests (see *power* above and Evidence).
5. **Prior record.** Done. "Every child lies in its parent's J = 0 component" moved from New to E-234 / E-239 / E-244 / E-254 / E-255 / E-268.

No item was left undone. No ball was re-run (the response changes none of their outputs; the referee reproduced all three).

## Pre-registration (written into this file before any count was run)

- P1. E-306's three balls hold exactly 2 key-kept perI J != 0 steps: 1 in 3060000 depth 4 (= E-233), 1 in 3304000 depth 5 (= E-239), 0 in 2400230 depth 3. **HIT** (1, 1, 0).
- P2. q_S equal on both; both already joined, so the outcome was predetermined [response: read "predicted by H-023 (OPEN) and E-306's record"]. **HIT.**
- P3. q_S equal on all 18 recorded key-kept J != 0 steps of the group (E-233, E-234 x10, E-239, E-244 x6), including E-244 8427 (never joined), confidence about 0.9. **HIT** 18 of 18; but 8427 is not open either: its child is isomorphic to 8063's child, which E-244 joined (E-255). The list missed E-254 (referee); with it, 19 of 19, the 19th not pre-registered.
- Predicted verdict: not zero, but no step with an unknown outcome. Holds.

## Claim

(a) Count (a rediscovery: E-306 Claim (c) and `toolsmith_qsball.txt` already state it; candidate totals are E-229/E-239). Rebuilt independently of E-306's loop, with every perI candidate tallied by outcome and `child()` failures counted rather than silently dropped: 3060000 depth 4 has 99 perI candidates (98 key-moved, 1 kept = E-233's step, path F1 R9 R8 R8, R8); 3304000 depth 5 has 423 (422 moved, 1 kept = E-239's index 2026); 2400230 depth 3 has none. 0 build failures, 0 illegal/vertex/err outcomes. So the count is 2, not 0.
(b) q_S on those steps and on all 17 other recorded key-kept J != 0 steps of the F-010 group: every recorded path rebuilds with `moves` moves only, every step is still perI-flagged and key-kept, and parent and child have the same (signature, |det S|, q_S) in 19 of 19 (11 in the 3060000 half, 8 in the 3304000 half incl. E-254 from 2400230). No step crosses halves; no child has an invariant of neither half.
(c) Why this is not a test of H-015. All 19 children were already joined to their parents' J = 0 components (E-233, E-234, E-239, E-244, E-254; 8427 via isomorphism to 8063, E-255; collected in E-254 and E-268). If H-023 holds (J = 0 steps are derived equivalences; OPEN), q_S is constant along those join chains, as it was on E-306's 14 815 ball steps. So q_S equality on these 19 was predicted by H-023 (OPEN) and E-306's record, and a q_S NO on a key-kept J != 0 step could only come from a child outside its parent's J = 0 component, or from a failure of H-023 on a join chain. No recorded child is outside. H-015 is testable at n = 9 by this instrument only in principle: the steps exist, but no recorded step has an outcome H-023 does not predict.
Not claimed: anything about H-015's truth (equal q_S does not imply equivalence); that q_S *can* say NO on a parent/child step (shown only for the halves, E-306); any shape common to the 19 steps; multi-vertex J != 0 steps; 3060000 beyond depth 5. Refuted by: a key-kept J != 0 step in these balls missing from my tally, or one of the 19 with different q_S on another machine.

## Evidence

| ball / source | parents scanned | perI candidates | key-kept | q_S equal | crosses halves |
|---|---|---|---|---|---|
| 3060000, depth <= 4 (E-306 ball) | 1584 | 99 (d3: 8, d4: 91) | 1 (E-233) | 1/1 | 0 |
| 3304000, depth <= 5 (E-306 ball) | 4046 | 423 (d3: 6, d4: 63, d5: 354) | 1 (E-239 2026) | 1/1 | 0 |
| 2400230, depth <= 3 (E-306 ball) | 130 | 0 | 0 | - | - |
| 3060000 depth-5 parents (E-234, replayed) | 10 paths | - | 10 | 10/10 | 0 |
| 3304000 depth-6 parents (E-244, replayed) | 6 paths | - | 6 | 6/6 | 0 |
| 2400230 depth-8 parent (E-254, replayed) | 1 path | - | 1 | 1/1 | 0 |
| total distinct recorded steps | | | 19 | 19/19 | 0 |

Shape key (printed by `replay`): step = the path's last move (same kind, vertex) 17 of 19; step not in the path 2 of 19 (E-234 4722, E-244 11123). Planted cross comparison (child vs other half's root invariant): NO 19 of 19, which is E-306's half-separation bit repeated, not a power control for the step test.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/121/skeptic_jcount.py replay             # 19 s (19 steps)
timeout 10m .venv/bin/python -W ignore workshop/rounds/121/skeptic_jcount.py ball 3060000 4    # 58 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/121/skeptic_jcount.py ball 3304000 5    # 212 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/121/skeptic_jcount.py ball 2400230 3    # 6 s
```
The script reuses `moves` (round 049), `perI` (round 033 via 049), `child` (round 082) and `qinv` (round 117) unchanged, so it shares their code with E-306; what is independent is the loop, the outcome tally, the path replay and the shape key.

## Prior record

E-306 (the count 2 and q_S on both: rediscovery, re-run and confirmed), E-229 and E-239 (candidate totals 99; 69 + 354), E-233, E-234, E-239, E-244, E-254 (the 19 steps and their joins). "Every n = 9 key-kept J != 0 child lies in its parent's J = 0 component": already recorded per step in E-234, E-239, E-244, E-255 (8427), collected in E-254's Result and restated in E-268 -- not new. H-023 (OPEN) is what predicts q_S equality along the joins. New, and small: q_S computed on the 17 deeper steps (predicted by H-023), the per-outcome tally showing E-306's loop dropped nothing, and the printed shape key (17 + 2) correcting E-239's "shared diamond shape" as a description of all steps.

## Code changed

None in the library. New script `workshop/rounds/121/skeptic_jcount.py` (response: E-254 added to STEPS, shape key printed); no tests touched.

## Next

- Chair / toolsmith: do not propose the 3060000 depth-6 ball as an overnight *test* of T10. It is also not overnight-sized: E-244 did the 3304000 depth-6 scan in commands of <= 10 min. If run, run it as a hit search with a q_S check on every key-kept child, pre-registering YES on all; only a child outside its parent's J = 0 component, or with a q_S of the other half, would be informative.
- The NO-capable route at n = 9 needs a step that leaves the J = 0 component: search key-kept J != 0 steps that perI does not flag (multi-vertex), or use the n = 12 47-group (H-032, E-314) where SNF(S) is a second instrument.
- Theorist: is there a reason every perI one-vertex key-kept J != 0 child returns into its parent's J = 0 component? It does not come from repeating a move at the same vertex: 4722 and 11123 do not.
