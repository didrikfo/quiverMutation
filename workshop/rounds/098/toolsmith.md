# Eight key-kept J != 0 edges (E-254, E-244, E-239), taken as undirected moves, give no meet of 2400230 (F/R depth 8) and 3304000 (depth 6) at radius 1 around the children (n = 9, one pair, capped radii)

author: toolsmith · round: 098 · kind: negative
thread: T10 · bears on: E-226, E-238, E-254, E-255, the J = 0 premise (T10/T5)
scope: only these 8 edges (as undirected moves), F/R radii 8 and 6, plus radius-1 balls of the 8 children; not the infinite move graph. n = 9; pair 2400230 (F/R ball depth 8, 20 210 nodes) vs 3304000 (depth 6, 15 401); J != 0 edges = the 8 key-kept perI steps found (1 from the 2400230 depth-8 parents, 7 from 3304000 parents of depth 5-6); each child gets its radius-1 J = 0 ball; lnaMoves checked at component level. 2400230 parents of depth < 8 are covered only by E-254's scan (none kept); 3304000 depth-6 children of depth-6 parents' own J != 0 steps beyond radius 1 not walked.

## Response to referee
Verdict minor revision; all four required items done. Outputs of jedges, meet, control were regenerated from the same pickles and committed as `toolsmith_j98_out.txt` (1.6 KB); they match the review's re-run (8 edges all LNA None; extras A 14 / B 94; sizes 20 224 / 15 495; 0 WL pairs, 0 iso tests; control 6 for 8063, 8427, 11123, d5, None for 8170, 8196, 4606 and the A child).
1. Output file committed: `workshop/rounds/098/toolsmith_j98_out.txt`. Done.
2. Undirected-edge assumption and the "only the A child's 14 nodes could hit" point are now in the Claim. Done.
3. Title and Scope narrowed to the 8 edges, radii 8 and 6, radius-1 balls of the children. "Does not close the gap" is dropped. Done.
4. Evidence table now says the A depth < 8 and B depth-6 scans are E-254's and E-244's and were not rerun (only the 6 kept parents and the depth-5 scan were replayed). Done.
Also taken from the review: the result follows from E-254/E-255 plus component membership (B children are iso to B nodes, the A child lies in A's component); non-LNA is real (9-10 arrows on 9 nodes, a vertex of degree 2), not a swallowed exception. The ball builds were not independently rerun by the referee; they stay the author's (sizes agree with E-244/E-254/E-255).

## Stated beforehand (meaning of the outcome)
- MEET (iso, sign-aware, any total): the J != 0 edges are a candidate bridge. It would show that the gap of E-238 is a missing J != 0 move, so J != 0 steps can lie in the F-010 equivalence of a certified pair. That supports "the J = 0 premise is too narrow as a move set" (J != 0 steps are derived steps there); it still is not a NO or a YES for the premise, because 2400230 ~ 3304000 is certified independently (E-238).
- NO MEET: a statement about this move set at these radii. It leaves the J = 0 premise untouched (both for and against), as E-250/E-255 already said, and it says the 8 recorded J != 0 edges are not the bridge.

## Claim
Assumption: every J != 0 edge is used as an undirected move (script docstring; not derived). With the 8 key-kept J != 0 edges so used (children plus radius-1 balls, +14 nodes on the A side, +94 on the B side) on top of the F/R balls, the isomorphism meet of 2400230(8+) and 3304000(6+) is empty: 0 Weisfeiler-Lehman bucket pairs, 0 isomorphism tests. The test is nearly vacuous by construction: A8 vs B6 already had 0 WL pairs (E-255), the 108 extra nodes are radius-1 neighbours of ball nodes, and the control shows the B children are iso only to B-ball nodes (depth 6), so only the A child's 14 nodes (depth 9-10 of A against a depth-6 B) had any chance of a hit. The result follows from E-254/E-255 plus component membership (B children lie in 3304000's J = 0 component, the A child in 2400230's). No extra node is an LNA (all 8 children have LNA name None), so lnaMoves adds nothing; the roots are in different lnaMoves components (E-255). Not claimed: that no J != 0 bridge exists (larger 2400230 balls are not built), that the gap of E-226/E-238 is closed or open as a property of the infinite move graph, nor anything about the J = 0 premise.

## Evidence
| step | result |
|---|---|
| coverage | A-side "none kept at depth < 8" is E-254's scan and B-side depth-6 scan is E-244's perI flag; neither rerun here (only the 6 kept parents replayed, plus the depth-5 scan of 2 958 parents) |
| rebuilt balls | A8 20 210 nodes, B6 15 401 (agree with E-254, E-244) |
| 2400230 side J != 0 | replay of E-254's path: parent at depth 8, one kept step (F at 3); child not LNA |
| 3304000 side J != 0 | the six E-244 parents (8063, 8170, 8196, 8427, 11123, 4606) each reproduce one kept step (R9, R7, R7, R7, R2, R7); depth-5 scan of all 2 958 parents gives 1 (R7, E-239's): 7 in all, none LNA |
| control (children vs own plain ball) | 4 of 7 B children are iso to a node of the plain depth-6 B ball at depth 6 (8063, 8427, 11123, depth-5 one), the other 3 (8170, 8196, 4606) lie beyond it and E-244 joined them at 7-8: matches E-244 (5 of 6 joined, total 6-8; 8427's join at 6 is E-255's iso with 8063's child). The A child is beyond A8's ball (E-254: total 9) |
| meet A(8)+14 vs B(6)+94 | 20 224 vs 15 495 nodes, WL-bucket pairs 0, iso tests 0, meet none (20 s) |
| lnaMoves | roots in different undirected components (as E-255); 0 extra LNA nodes |

The empty bucket intersection is the same signature the E-226/E-254/E-255 meets show: the two balls hold no graphs with equal WL hash at all, so a larger iso budget could not change this verdict. Caveat: the meet is a necessary-test meet; with 0 pairs nothing is "decided" beyond the move set.

## Reproduction
Builds are not committed (pickles rebuilt in the scratchpad; two parallel runs of each took 8.4 + 5 min):
```
DEPTH=8 timeout 10m .venv/bin/python workshop/rounds/087/toolsmith_d6scan.py build 2400230 A8.pkl --budget-hours 0.14   # exit 2 -> rerun, 2 commands
DEPTH=6 timeout 10m .venv/bin/python workshop/rounds/087/toolsmith_d6scan.py build 3304000 B6.pkl --budget-hours 0.14   # 2 commands
timeout 10m .venv/bin/python workshop/rounds/098/toolsmith_j98.py jedges A8.pkl B6.pkl J.pkl     # 84 s
timeout 10m .venv/bin/python workshop/rounds/098/toolsmith_j98.py meet A8.pkl B6.pkl J.pkl       # 26 s
timeout 10m .venv/bin/python workshop/rounds/098/toolsmith_j98.py control A8.pkl B6.pkl J.pkl
```
Outputs of the last three: `workshop/rounds/098/toolsmith_j98_out.txt` (committed).
```
```

## Prior record
Largely known: E-254 (A-side child no meet with 3304000 radius 5; 2400230(6) vs 3304000(5)), E-255 (F/R + one lnaMoves layer at 8 vs 6, 0 common components, 8427's child iso to 8063's), E-238 (certified equivalent). New here only: all eight J != 0 edges together on the same union at 8 vs 6, with radius-1 balls, and the observation that every such child is non-LNA, so a bridge would have to come from F/R balls of J != 0 children, not lnaMoves. Component-level argument: the recorded B-side children join 3304000's own J = 0 component and the A-side child joins 2400230's, so in the infinite union graph these edges connect nothing new; only the finite-radius ball statement above is measured.

## Code changed
New `workshop/rounds/098/toolsmith_j98.py` (modes --plan, jedges, meet, control). No library change, no tests needed.

## Next
Not an overnight case: more radius is the same F/R-ball statement and every recorded J != 0 edge sits inside a single move component, so it cannot join them. The only edge that could would be a J != 0 step of a node at 2400230 depth 9-10 (E-254's 55 000 / 150 000-node proposal, hours, unlikely: 1 kept step per 20 210 parents). Skeptic/theorist: a certificate for a non-LNA pair that is not a move set is what T10 waits for; I would close the 2400230 calibration here.
