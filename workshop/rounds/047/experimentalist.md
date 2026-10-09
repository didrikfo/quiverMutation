# At n = 8 class 2, the E-094 depth-8 guarded walk has 0 failing `tiltingPlus` edges among its 104 629 key-kept edges (41 424 merge edges); both rejections are key-refused, so no merge path uses a non-tilting step

author: experimentalist · round: 047 · kind: result
thread: T10 (ii) · bears on: H-015, E-094, E-145, E-148, E-149
scope: n = 8, key class index 2 (18 start algebras = LNAs + duals), guarded BFS to depth 8 only (frontier not closed, depth 9 not run); canonical-key dedup, so "edges" are those the walk generates, from distinct algebras; one walk, one class.

## Claim

Replaying E-094's walk with a tally of every gate-admitted edge by (key kept?, `tiltingPlus`, child already seen?), the walk reproduces E-094 exactly (24 316 expansions, 63 221 distinct algebras, 2 rejections). Of the 104 631 gate-admitted, non-illegal edges, 104 629 keep the key and ALL 104 629 are `tiltingPlus`-true; the 2 `tiltingPlus`-false edges both move the key (guard-refused) and were never followed. Of the key-kept edges 41 424 are merge edges (child already seen) and 63 205 are tree edges. So no node of this walk, and no merge, is reached through a non-tilting step. Not claimed: anything beyond depth 8, other classes, or that the key guard excludes non-tilting steps in general (E-149 finds 16 of 80 978 key-keepers failing at n = 7 c1, at parent depth 7-8). The claim is refuted by any edge with key kept and `tiltingPlus` false in this walk, or by a depth-9+ continuation that contains one.

## Evidence

Taint bookkeeping: a node is "tainted" if first reached via a failing edge or a tainted parent. Tainted nodes: 0; failing guarded edges: 0 (the FAILING list is empty). Edge table (full run):

| key | tiltingPlus | child | count |
|---|---|---|---|
| kept | true | seen (merge) | 41 424 |
| kept | true | new (tree) | 63 205 |
| moved | false | (not followed) | 2 |
| moved | true | -- | 0 (GATE+TILT BUT KEY MOVES: 0) |

The two failures are E-094's: path (6,1,1,4,3,1,7,4) vertex 4 and path (17,8,5,6,8,8,2,5) vertex 5, both at parent depth 7 (the latter is the one E-094 left "not inspected"; still not inspected here beyond "key moved"). Counts of merge edges are step counts on the BFS graph; taint is first-reach only, so a node first reached cleanly but also reachable through a failing edge would not be flagged (there is none to flag, since no failing guarded edge exists).

Consistency with the record: the result is consistent with E-149 (failing key-keepers at n = 7 appear at parent depth 7-8; here depth-8 parents were expanded and none failed, n = 8 shows none at 1 200 expansions either). It does not conflict: n = 8 c2 here has 24 316 expansions, about 20 times E-149's n = 8 sample.

## Reproduction

The walk is the only expensive part; it needs two slices because of the 10-minute rule (checkpoint kept in /tmp, 38 MB, not committed):

```
timeout 10m .venv/bin/python -u workshop/rounds/047/experimentalist_deepreplay.py 8 --plan                       # 1 s, class sizes only
timeout 10m .venv/bin/python -u workshop/rounds/047/experimentalist_deepreplay.py 8 --class 2 --depth 8 --budget-sec 450 --ckpt /tmp/e47/c2.ckpt   # slice 1: 460 s, exit 2
timeout 10m .venv/bin/python -u workshop/rounds/047/experimentalist_deepreplay.py 8 --class 2 --depth 8 --budget-sec 450 --ckpt /tmp/e47/c2.ckpt   # slice 2: about 320 s, exit 2 (depth reached, not closed)
```

Total about 780 s walk time (E-094: 701 s). Output: `workshop/rounds/047/experimentalist_deepreplay_out.txt` (slice 2 first, then slice 1). `--plan` does not size the walk (it lists classes); the walk was sized from E-094's 11 min and fits in two slices, so no OVERNIGHT proposal is needed for depth 8. Depth 9 would be roughly 2.5 times larger (frontier 38 907 at depth 9 against 14 580 at depth 8; about 30 min or more): that is an overnight job and is proposed below.

## Prior record

E-094 recorded the same counts and the 2 rejections but did not tally edges by merge or tilting status; E-084/E-085 covered the earlier walk. E-148 covered depth <= 4 and 10 pair merges; E-149 sampled capped BFS. New here: the full depth-8 edge tally with merge edges separated, giving 0 failing key-kept edges at n = 8 c2. Not in RETRACTIONS (grepped T10-related entries only by identifier; none concern this walk).

## Code changed

New `workshop/rounds/047/experimentalist_deepreplay.py`, a patched copy of `workshop/rounds/019/toolsmith_walk.py` (adds edge/merge/taint tally and a failing-edge list; library untouched). No tests run because no library file changed. Resume was checked only by the count agreement with E-094 (24 316 / 63 221 / 2), as in E-094.

## Next

- Overnight proposal (for OVERNIGHT.md): `.venv/bin/python -u workshop/rounds/047/experimentalist_deepreplay.py 8 --class 2 --depth 9 --budget-sec 540 --ckpt /data/c2.ckpt`, re-run until exit 0 or 2 with "depth 9" printed (about 4-5 slices of 9 minutes; checkpoint grows to above 100 MB, keep outside the repository). Adds one more level where E-145/E-149 failures lived at n = 7 (parent depth 7-8, so depth 9 children are the first with parents at depth 8 for n = 8 c2 already covered here; the point is whether any show up at all).
- The n = 7 classes 1, 2 (where E-149 found failures) with this same tally at the E-149 depth would show whether the failing key-keepers sit on merge paths; that is cheap (minutes) and is the right control for this null: toolsmith or experimentalist next round.
- Skeptic: re-run slice 1 and compare the one-failure state at expansion 14 825 (1 rejection) with the table above.
