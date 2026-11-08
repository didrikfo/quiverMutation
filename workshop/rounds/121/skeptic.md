# At n = 9 (F-010 key group, balls of depth <= 5 / 6 only) every recorded key-kept J != 0 step keeps q_S: E-306's balls hold 2 such steps (1 + 1 + 0, as E-306 says), the record holds 18, and all 18 children were already known to sit in their parent's J = 0 component, so q_S cannot say NO on any of them; H-015 has no informative step at n = 9 yet

author: skeptic · round: 121 · kind: negative (mostly rediscovery of E-306's count, plus the q_S test on all 18 recorded steps)
thread: T10 · bears on: H-015, H-032, H-023
scope: n = 9, F-010 key group only; E-212 `moves` balls: 3060000 depth 4 (parents depth <= 4), 3304000 depth 5 (parents depth <= 5), 2400230 depth 3, plus the 16 further recorded steps replayed from their paths (E-234: 3060000 depth-5 parents; E-244: 3304000 depth-6 parents); only perI-flagged one-vertex F/R J != 0 candidates (multi-vertex J != 0 not flagged by perI are not seen); 3060000 depth-5 ball not rebuilt (its 10 steps replayed from paths), 3060000 depth 6 not run.
power: planted cross steps (each of the 18 children compared with the OTHER half's root invariant, the halves certified inequivalent by F-010/E-238): q_S says NO 18 of 18. On the real steps the test could only fire on a child outside its parent's J = 0 component, and none of the 18 is (see Claim), so its power on H-015 here is nil.

## Pre-registration (written into this file before any count was run)

- P1. E-306's three balls hold exactly 2 key-kept perI J != 0 steps: 1 in 3060000 depth 4 (= E-233), 1 in 3304000 depth 5 (= E-239), 0 in 2400230 depth 3. **HIT** (1, 1, 0).
- P2. q_S equal on both; both already joined, so the outcome was predetermined. **HIT.**
- P3. q_S equal on all 18 recorded key-kept J != 0 steps of the group (E-233, E-234 x10, E-239, E-244 x6), including E-244 8427 (never joined), confidence about 0.9. **HIT** 18 of 18; but 8427 is not open either: its child is isomorphic to 8063's child, which E-244 joined (E-255; I had missed this before registering).
- Predicted verdict: not zero, but no step with an unknown outcome. Holds.

## Claim

(a) Count (a rediscovery: E-306 Claim (c) and `toolsmith_qsball.txt` already state it). Rebuilt independently of E-306's loop, with every perI candidate tallied by outcome and `child()` failures counted rather than silently dropped: 3060000 depth 4 has 99 perI candidates (98 key-moved, 1 kept = E-233's step, path F1 R9 R8 R8, R8); 3304000 depth 5 has 423 (422 moved, 1 kept = E-239's index 2026); 2400230 depth 3 has none at all. 0 build failures, 0 illegal/vertex/err outcomes, candidate totals equal E-229/E-239 (99; 69 + 354). So the count is 2, not 0.
(b) q_S on exactly those steps, and on all 16 other recorded key-kept J != 0 steps of the F-010 group: every recorded path rebuilds with `moves` moves only, every step is still perI-flagged and key-kept, and parent and child have the same (signature, |det S|, q_S) in 18 of 18 (11 in the 3060000 half, 7 in the 3304000 half). No step crosses halves; no child has an invariant of neither half.
(c) Why this is not a test of H-015. All 18 children were already joined to their roots by J = 0 chains (E-233, E-234, E-239, E-244; 8427 via isomorphism to 8063, E-255), and q_S has been constant on all 14 815 J = 0 steps computed (E-306). So q_S equality on these 18 was implied before it was computed; a q_S NO on a key-kept J != 0 step can occur only if the child lies outside its parent's J = 0 component. Every key-kept J != 0 child found at n = 9 so far (18 of 18, all the "R at v twice" diamond shape of E-239) lies inside it. H-015 is therefore testable at n = 9 by this instrument only in principle: the instrument can say NO (planted 18/18), the steps exist, but no known step has an open outcome.
Not claimed: anything about H-015's truth (a YES from q_S has no power: equal q_S does not imply equivalence); multi-vertex J != 0 steps; 3060000 beyond depth 5. Refuted by: a key-kept J != 0 step in these balls missing from my tally, or one of the 18 with different q_S on another machine.

## Evidence

| ball / source | parents scanned | perI candidates | key-kept | q_S equal | crosses halves |
|---|---|---|---|---|---|
| 3060000, depth <= 4 (E-306 ball) | 1584 | 99 (d3: 8, d4: 91) | 1 (E-233) | 1/1 | 0 |
| 3304000, depth <= 5 (E-306 ball) | 4046 | 423 (d3: 6, d4: 63, d5: 354) | 1 (E-239 2026) | 1/1 | 0 |
| 2400230, depth <= 3 (E-306 ball) | 130 | 0 | 0 | - | - |
| 3060000 depth-5 parents (E-234, replayed) | 10 paths | - | 10 | 10/10 | 0 |
| 3304000 depth-6 parents (E-244, replayed) | 6 paths | - | 6 | 6/6 | 0 |
| total distinct recorded steps | | | 18 | 18/18 | 0 |

Planted cross steps (child vs other half's root invariant): NO 18 of 18.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/121/skeptic_jcount.py replay             # 18 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/121/skeptic_jcount.py ball 3060000 4    # 58 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/121/skeptic_jcount.py ball 3304000 5    # 212 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/121/skeptic_jcount.py ball 2400230 3    # 6 s
```
The script reuses `moves` (round 049), `perI` (round 033 via 049), `child` (round 082) and `qinv` (round 117) unchanged, so it shares their code with E-306; what is independent is the loop, the outcome tally, and the path replay.

## Prior record

E-306 (the count 2 and q_S on both: this part is a rediscovery, re-run and confirmed), E-229 (candidate totals), E-233, E-234, E-239, E-244 (the 18 steps and their joins), E-255 (8427 isomorphic to 8063). New: q_S on the 16 deeper steps (not recorded before), the per-outcome tally showing no silent drop, and the observation that every n = 9 key-kept J != 0 child lies in its parent's J = 0 component, which is what makes the q_S test uninformative on them.

## Code changed

None in the library. New script `workshop/rounds/121/skeptic_jcount.py`; no tests touched.

## Next

- Chair / toolsmith: do not propose the 3060000 depth-6 ball as an overnight *test* of T10. It is also not overnight-sized: E-244 did the 3304000 depth-6 scan in commands of <= 10 min. If run, run it as a hit search with a cheap q_S check on every key-kept child, and pre-register the expected outcome (YES on all, diamond shape); only a child outside its parent's J = 0 component, or with a q_S of the other half, would be informative.
- The NO-capable route at n = 9 needs a step that leaves the J = 0 component: search for key-kept J != 0 steps not of the "R at v twice" diamond shape (e.g. multi-vertex J != 0 steps, which perI does not flag), or use the n = 12 47-group (H-032, E-314) where SNF(S) is a second instrument.
- Theorist: is it a theorem that a perI one-vertex step repeated at the same vertex (the diamond shape) returns into the parent's J = 0 component? If so, the n = 9 hits are uninformative by construction and the count of 18 is a count of a known shape.
