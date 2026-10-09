# At n = 10, merges.py (depths 3, 4 all 122 members; depth 5 only 19) finds no link, and the one recorded merge replayed (34504030 -> 50505000, 5 paths of 7) uses no J != 0 step

author: experimentalist · round: 051 · kind: negative
thread: T10 (ii) · bears on: F-037, F-041, E-032, E-033, E-145, E-148, E-154
scope: n = 10 only; the 122 members of the 12 leftover orbits in groups with >= 2 orbits (`merges.py` default, no `--all-groups`); depth 3 and 4 complete, depth 5 only 19 of 122 (time); links found by the run: 0. Of the recorded n = 10 merges, only F-037 (34504030 -> 50505000) was replayed; the second (05040330 -> 33460000) was NOT.

## Claim

(1) `merges.py 10 --depths 3 4 5 --witness` (4 jobs, 10-minute chunks, resumed from one checkpoint) produced no link at all: depth 3 and 4 over all 122 searched members, depth 5 over 19. So the run contains no witness link and cannot show a J != 0 dependence; it is consistent with E-032, where both n = 10 merges appear only at one-sided depth 6 (group A) and 7 (F-037), beyond depth 5. The assignment's "depth 5" therefore cannot answer the question on its own. This is a cap, not a verdict: depth 5 is 84 percent undone.
(2) A real witness link was replayed instead: the five F-037 paths 34504030 -> 50505000 (vertex sequences from F-037; 7 mutations each) were replayed edge by edge from 34504030 (Kupisch row [3,4,5,0,4,0,3,0]): all 35 edges are gate-admitted, J = 0 under `tiltingPlus`, and key-keeping; each ends on the LNA row [5,0,5,0,5,0,0,0] (via `lnaMoves.asRelLengths`, i.e. up to relabelling, as `merges.py` decides it). So this recorded n = 10 merge does not depend on a J != 0 step along these 5 paths (5 of 19 recorded).
Not claimed: that no n = 10 merge uses a J != 0 step. The group-A merge (05040330 -> 33460000, 7 steps, E-033) is unreplayed here; E-033 checked its steps for gate, illegal relation, parallel arrow and key but not `tiltingPlus`. A refuting result: a replay of an A path with a `tiltingPlus` failure, or a merge found at depth <= 7 only through such an edge.

## Evidence

| depth | searches done | distinct members | with link | alarms | CPU-s |
|---|---|---|---|---|---|
| 3 | 122 | 122 | 0 | 0 | 1271 |
| 4 | 122 | 122 | 0 | 0 | 3575 |
| 5 | 19 | 19 | 0 | 0 | 1772 |

Cost: depth 3 about 12 s, depth 4 about 30 s, depth 5 about 80-150 s per search (sizing from this run; depth 5 over all 122 is about 4 CPU-hours, about 1 h wall on 4 cores; depth 6 about 5x, depth 7 for A and F is the hours-scale of E-032, slowest single search 9691 s at depth 8). The default call plans nothing; `--summary` shows the groups: group A (69 + 42 members), the 4-orbit group (4,2,2,1), F (1 + 1), and singleton groups skipped without `--all-groups`.
Side test, `experimentalist_a_meet.py`: depth-3 reach (with duals) of all 111 group-A members, hash-joined across the two orbits on labelled quiver keys: 0 meetings of total 6, as expected for a 7-step link (a 3 + 3 meeting needs total <= 6 and labels equal). It is a null that agrees with E-033, not a new result. A depth-4 meet from 05040330 and 33460000 on labelled keys also found nothing (that script was dropped): the link's target is an LNA only up to relabelling, so labelled meetingPoints may miss it even at total 7 (untested guess); recognising the end by `asRelLengths` or an isomorphism key is the fix.
F-037 replay: 5 paths, 7 edges each, gate 35, J = 0 35, key 35 of 35, all end on the target row.

## Reproduction

```
timeout 10m .venv/bin/python -u merges.py 10 --depths 3 4 5 --witness --jobs 4 --budget-hours 0.13 --checkpoint workshop/rounds/051/experimentalist_merges10.jsonl   # repeat 4 times (resumes); about 10 min each
.venv/bin/python workshop/rounds/051/experimentalist_merges10_summary.py     # table above; seconds
timeout 10m .venv/bin/python workshop/rounds/051/experimentalist_f037replay.py   # about 10 s
timeout 10m .venv/bin/python workshop/rounds/051/experimentalist_a_meet.py       # about 5 min
```
The checkpoint (94 KB) and logs are in `workshop/rounds/051/experimentalist_merges10*.{jsonl,log}`.

## Prior record

E-032/E-033/F-037: the two n = 10 merges and the depths (6-7) at which merges.py finds them; E-154: `--witness` untested on a real link (still untested: no link appeared). E-148/E-154: n = 8 F-041 merges J = 0. New here: n = 10 depth <= 4 complete null and depth 5 partial in this configuration; `tiltingPlus` on all 35 edges of F-037's paths. Mild point against reading "depth 5" as sufficient: the recorded links are deeper.

## Code changed

None in the library or `merges.py`. New scripts in `workshop/rounds/051/`: `experimentalist_merges10_summary.py`, `experimentalist_f037replay.py`, `experimentalist_a_meet.py`. No tests run (no code under test changed).

## Next

Proposal for `OVERNIGHT.md`: `merges.py 10 --depths 5 6 7 --witness --jobs 7 --checkpoint workshop/rounds/051/experimentalist_merges10.jsonl` (resumes; depth 5 about 1 h, depth 6 several hours, depth 7 for A and F at E-032 cost), then `tiltingPlus` replay of each recorded witness (the `witnesses` field holds `{start, path}`; start 1 = dual). Cheaper: one search `05040330` depth 7 with witness, then replay. Toolsmith: a labelled-vs-isomorphism meet in the middle for LNA targets (`meetingPoints` compares labelled keys).
