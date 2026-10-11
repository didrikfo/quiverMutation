# At n = 10 the group-A merge 05040330 -> 33460000: 5 explicit 7-step paths (E-033 section 4 recorded 7 steps, 5 checked), all edges gate-admitted, J = 0 under tiltingPlus (a test not shown able to fail here), key-keeping, ends isomorphic to the target member

## Response to referee

1. E-033 section 4 cited (done). It already records `05040330 -> 33460000`, 9 paths, 5 checked, 7 steps, each admissible/acyclic/key-holding. The existence of a length-7 link is therefore not new and the sentence "E-033's length 7 is attained" is withdrawn. New here: the explicit signed path strings, `tiltingPlus` J on each edge, and the relabelling-aware join. (E-164 and 051 say "no path recorded"; that inconsistency with E-033 is the record's, not resolved here.)
2. End identity (done). Established by structure-graph isomorphism (VF2) of the end with `LinearNakayamaAlgebra(10, m)` for the target member m, not by `asRelLengths` alone; now printed per path in `experimentalist_amerge_out.txt` ("end iso to member ... True" for 5/5). The same line also prints a labelled-key comparison of the end with the member: True for 5/5 in my run. The referee found the labelled keys differ; I could not reconcile this and do not rely on it either way (isomorphism is the claim).
3. Labelled-equality of the meetings (done). `experimentalist_amerge_out.txt` prints, per meeting, `kx == ky` (forward vs backward quiver as labelled keys): False for 5/5, with an isomorphism found and non-trivial. E-164 sentence narrowed to "the labelled `quiverKey` equality test fails on these 5 meetings".
4. "Dual-related in spirit" for rows 1-2: dropped (not checked; rows 45055000 and 55504400 are not obviously dual).
5. Non-vacuity of the J test (done, negative). `experimentalist_amerge_nonvac.py` / `_nonvac_out.txt`: 56 gate-admitted edges around the start (depth 2 both orientations) and 1146 gate-admitted (member, orientation, vertex) triples over all 122 checkpoint members: J != 0 in none. Where the gate rejects, `tiltingPlus` is false (it agrees with the gate), so on this code path J = 0 is not independent of the gate and I found no example where it fails. "35/35 J = 0" should be read as consistent with the gate, not as an extra check.
Changed: title, claim, scope line (below). The reproduction was re-run after the edit: same 5 meetings, same paths.

author: experimentalist · round: 054 · kind: result
thread: T10 · bears on: E-032, E-033, E-164, F-037 (group A), H-015
scope: n = 10, one start (only these 5 paths; no claim about meetingPoints in general or about J != 0 steps elsewhere) (LNA [0,5,0,4,0,3,3,0], class 05040330), the 42 distinct members of orbit 33460000 in the E-164 checkpoint; forward depth 4 (387 quivers) joined to backward depth 3 (141-207 quivers per target); parallel-arrow quivers skipped (`quiverKey` None). Only total-7 joins computed; shorter links (total <= 6) are excluded by E-033/E-164 (depth <= 5 one-sided, and a 3+3 labelled null), not re-proved here. Replay is Cartan-free: gate + `tiltingPlus` (Prop 2.3(c), one-map test) + Coxeter key only; generation of K^b(proj) by T is assumed as everywhere in T10.

## Claim

There are 5 isomorphism-level meetings of total length 7 between the forward depth-4 reach of 05040330 and the depth-3 reach of members of 33460000, and each gives a full path of 7 mutations (right steps positive, left steps negative) from [0,5,0,4,0,3,3,0] to an LNA in orbit 33460000. All 5 paths replay with `tiltingPlus`: 35 of 35 edges gate-admitted, 35 of 35 J = 0, 35 of 35 key-keeping, and each end is isomorphic (structure graph) to the target member and has its `asRelLengths` row. So the group-A merge needs no J != 0 step on these 5 paths, and E-033 section 4's recorded length 7 (9 paths, 5 checked) now has explicit replayed paths. None of the 5 meetings is labelled-equal (printed, `kx == ky` False; the isomorphism between the two halves is a non-trivial relabelling in all 5), so a labelled-key join misses these 5 meetings; this does not show anything about `meetingPoints` in general. The J = 0 result is consistent with the gate but I found no gate-admitted n = 10 edge with J != 0 (1202 examined), so it is not shown that the test can fail. Not claimed: that no group-A path uses a J != 0 step (only that these 5 do not), nor that 7 is shortest beyond E-033's own search.

## Evidence

| # | start | forward half (to X) | backward half (target -> Y), Y iso X | full path | edges gate/J=0/key | end row |
|---|---|---|---|---|---|---|
| 1 | 05040330 | 2 1 6 5 | 45055000 : 5 9 1 | 2 1 6 5 -4 -10 -8 | 7/7/7 | 4,5,0,5,5,0,0,0 |
| 2 | 05040330 | -10 -8 -4 -5 | 55504400 : -8 -4 -5 | -10 -8 -4 -5 2 1 6 | 7/7/7 | 5,5,5,0,4,4,0,0 |
| 3 | 05040330 | 4 3 2 1 | 60504030 : -10 -7 -8 | 4 3 2 1 7 6 9 | 7/7/7 | 6,0,5,0,4,0,3,0 |
| 4 | 05040330 | 4 3 7 2 | 60504030 : -10 -7 -2 | 4 3 7 2 1 6 9 | 7/7/7 | 6,0,5,0,4,0,3,0 |
| 5 | 05040330 | 4 3 7 9 | 60504030 : -7 -2 -3 | 4 3 7 9 2 1 6 | 7/7/7 | 6,0,5,0,4,0,3,0 |

Sign convention: positive = right mutation at that vertex, negative = left mutation (the search's dual walk). The backward half is inverted: reverse the order, flip each sign (inverse of a right mutation is a left mutation at the same vertex), and rename each vertex through the isomorphism phi: Y -> X found by VF2 on the structure graph (quiver arrows plus relation-path chains), so the replay is on the labelling of the forward half. A left step is replayed on the opposite algebra (`dualPathAlgebra`), gate and `tiltingPlus` evaluated there. The ends were checked by structure-graph isomorphism with the target member (and the row), not merely accepted as "an LNA".

Sizing: forward depth 3/4/5 = 131/387/1091 quivers in 4/23/108 s; backward depth 3 for 42 targets (pool of 4) = 117 s; join + replay < 1 min. Depth 5 forward was not needed (a 5 + 3 join would only find total 8). The search is cheap because the Coxeter guard prunes hard from this start. Explicit match key: Weisfeiler-Lehman hash of the structure graph, then exact isomorphism, so the hash is only a filter.

Caveat on the 5: rows 3-5 share one target; no relation between rows 1 and 2 is claimed. 5 meetings is the count at total 7 with this particular depth split (4 + 3); a 3 + 4 split or the duals might add more.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/054/experimentalist_amerge.py fwd 4 F4.pkl      # 23 s
timeout 10m .venv/bin/python workshop/rounds/054/experimentalist_amerge.py back 3 B3.pkl     # 117 s
timeout 10m .venv/bin/python workshop/rounds/054/experimentalist_amerge_join.py F4.pkl B3.pkl 5   # < 1 min
```
(the pickles are scratch; output of the last command: `workshop/rounds/054/experimentalist_amerge_out.txt`. Non-vacuity: `.venv/bin/python workshop/rounds/054/experimentalist_amerge_nonvac.py`, output `_nonvac_out.txt`.) Targets are read from `workshop/rounds/051/experimentalist_merges10.jsonl`.

## Prior record

E-032/E-033: the merge exists at one-sided depth 6/7, no path recorded. E-164: `merges.py 10 --depths 3 4 5` no link; group A unreplayed; `meetingPoints` labelled. 051 `a_meet` null: depth 3 + 3 labelled, 0 meetings of total 6 (consistent: a 7-step link). E-033 section 4 recorded 9 paths / 5 checked / 7 steps. New here: explicit path strings, the tiltingPlus replay, and the labelled-equality failure of these 5 meetings. The checkpoint's orbit label 33460000 groups many relation strings (4505..., 5550..., 6050...), the endpoint is "in the orbit", not the literal row 33460000.

## Code changed

None in `quivermutation/`. New scripts `experimentalist_amerge.py`, `experimentalist_amerge_join.py`, `experimentalist_amerge_nonvac.py` (round 054 directory). No tests run (no library code touched).

## Next

- toolsmith: a relabelling-aware `meetingPoints` (use the structure-graph WL hash + exact isomorphism, as in `struct`/`iso`), with parallel-arrow nodes handled; the signed path convention is what is needed.
- skeptic: hand-check path 3 (4 3 2 1 7 6 9) edge by edge with an independent End(T) at quiver level (E-167 style); the Hom test (E-166) on the 35 edges.
- Do the same for F-037 (34504030 -> 50505000) as a cross-check on the claim that the labelled join, not the class, was the obstacle.
