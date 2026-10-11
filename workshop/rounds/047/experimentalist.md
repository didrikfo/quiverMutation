# At n = 8 class 2, no key-kept edge among the 104 629 generated in the E-096 depth-8 guarded walk fails `tiltingPlus` (the 2 failures are key-refused)


## Response to referee

Verdict was minor revision; three required items, all addressed. The title and claim are narrowed as the referee suggests: the claim is now only "no key-kept edge among the 104 629 generated in the depth-8 walk fails `tiltingPlus`". I no longer state the stronger "no merge path uses a non-tilting step" as a result; it is a corollary of the narrowed claim for this walk's edges, not a separate finding.

1. Tree edges + starts vs distinct algebras (off by 2; referee's slice off by 1). Explained from the code and the saved checkpoint (no new run; one 5 s read of /tmp/e47/c2.ckpt). In `add()`, a child whose `canonicalKey` is `None` (relabeling cap exceeded) is appended to the queue as "new" but is never put in `seen`, because a `None` key cannot be recognised (the library's safe direction). Check: the per-depth queue sizes (18, 78, 194, 444, 1030, 2344, 5628, 14580, 38907) sum to 63 223 = 18 starts + 63 205 tree edges, exactly; `len(seen)` is 63 221. Recomputing the keys of the 38 907 depth-8 queue entries in the checkpoint gives exactly 2 with key `None` and 0 keyed entries missing from `seen`. So the 2 are two keyless depth-8 children counted as tree edges. The 18 starts are all distinct and keyed (depth-1 expansions = 18). The referee's off-by-1 at 14 674 expansions is consistent (one of the two keyless children was generated before that point); I did not check which slice-1 child it is.
   What it affects: "63 221 distinct algebras" is really 63 221 distinct keyed algebras plus 2 unkeyed nodes of unknown identity (they could coincide with each other or with a keyed node under a different presentation; not tested). The tree figure 63 205 includes those 2, and the merge figure 41 424 is a lower bound: an edge into a keyless node can never be counted as a merge. So the split is "tree 63 205 (at most 2 too many) / merge 41 424 (possibly a few too few)". The total 104 629, the tilting tally and the claim do not depend on the split. The 63 221 / 24 316 / 2 agreement with E-096 is unaffected (same convention there). Not a discrepancy left open; only the identity of the 2 keyless nodes is unknown.

2. Taint bookkeeping demoted. "Tainted nodes: 0" is removed from the Evidence as a finding. With 0 failing guarded edges there is nothing to taint, so it adds nothing beyond the edge table. The script still prints the tally for completeness; it is not evidence.

3. Resume and independent check. The final tally spans a checkpoint resume (slice 1 ended at depth 8 pos 5089, expansion 14 825; slice 2 resumed from the pickled state). Checking of the resume is the count agreement with E-096 (24 316 / 63 221 / 2) only. The skeptic independently re-ran slice 1 only (stopped at expansion 14 674; same depth 1-7 lines, same one rejection, 0 failing, 0 key-moves-but-tilting); that matches the author's slice-1 partial tally. Slice 2 (about 320 s) was not independently re-run: the final 41 424 / 63 205 / 2 rest on the author's output file plus the E-096 count agreement.

author: experimentalist · round: 047 · kind: result
thread: T10 (ii) · bears on: H-015, E-096, E-147, E-150, E-151
scope: n = 8, key class index 2 (18 start algebras = LNAs + duals), guarded BFS to depth 8 only (frontier not closed, depth 9 not run); canonical-key dedup, so "edges" are those the walk generates, from distinct algebras; one walk, one class.

## Claim

Replaying E-096's walk with a tally of every gate-admitted edge by (key kept?, `tiltingPlus`, child already seen?), the walk reproduces E-096 exactly (24 316 expansions, 63 221 distinct algebras, 2 rejections). Of the 104 631 gate-admitted, non-illegal edges, 104 629 keep the key and ALL 104 629 are `tiltingPlus`-true; the 2 `tiltingPlus`-false edges both move the key (guard-refused) and were never followed. Of the key-kept edges 41 424 are merge edges (child already seen) and 63 205 are tree edges. So every key-kept edge this walk generated is tilting. Not claimed: anything beyond depth 8, other classes, or that the key guard excludes non-tilting steps in general (E-151 finds 16 of 80 978 key-keepers failing at n = 7 c1, at parent depth 7-8). The claim is refuted by any edge with key kept and `tiltingPlus` false in this walk, or by a depth-9+ continuation that contains one.

## Evidence

Failing guarded edges: 0 (the FAILING list is empty; the script's taint tally is vacuous for that reason and is not offered as evidence). Edge table (full run; the merge/tree split is approximate by 2 keyless children, see the response above):

| key | tiltingPlus | child | count |
|---|---|---|---|
| kept | true | seen (merge) | 41 424 |
| kept | true | new (tree) | 63 205 |
| moved | false | (not followed) | 2 |
| moved | true | -- | 0 (GATE+TILT BUT KEY MOVES: 0) |

The two failures are E-096's: path (6,1,1,4,3,1,7,4) vertex 4 and path (17,8,5,6,8,8,2,5) vertex 5, both at parent depth 7 (the latter is the one E-096 left "not inspected"; still not inspected here beyond "key moved"). Counts of merge edges are step counts on the BFS graph.

Consistency with the record: the result is consistent with E-151 (failing key-keepers at n = 7 appear at parent depth 7-8; here depth-8 parents were expanded and none failed, n = 8 shows none at 1 200 expansions either). It does not conflict: n = 8 c2 here has 24 316 expansions, about 20 times E-151's n = 8 sample.

## Reproduction

The walk is the only expensive part; it needs two slices because of the 10-minute rule (checkpoint kept in /tmp, 38 MB, not committed):

```
timeout 10m .venv/bin/python -u workshop/rounds/047/experimentalist_deepreplay.py 8 --plan                       # 1 s, class sizes only
timeout 10m .venv/bin/python -u workshop/rounds/047/experimentalist_deepreplay.py 8 --class 2 --depth 8 --budget-sec 450 --ckpt /tmp/e47/c2.ckpt   # slice 1: 460 s, exit 2
timeout 10m .venv/bin/python -u workshop/rounds/047/experimentalist_deepreplay.py 8 --class 2 --depth 8 --budget-sec 450 --ckpt /tmp/e47/c2.ckpt   # slice 2: about 320 s, exit 2 (depth reached, not closed)
```

Total about 780 s walk time (E-096: 701 s). Output: `workshop/rounds/047/experimentalist_deepreplay_out.txt` (slice 2 first, then slice 1). `--plan` does not size the walk (it lists classes); the walk was sized from E-096's 11 min and fits in two slices, so no OVERNIGHT proposal is needed for depth 8. Depth 9 would be roughly 2.5 times larger (frontier 38 907 at depth 9 against 14 580 at depth 8; about 30 min or more): that is an overnight job and is proposed below.

## Prior record

E-096 recorded the same counts and the 2 rejections but did not tally edges by merge or tilting status; E-086/E-087 covered the earlier walk. E-150 covered depth <= 4 and 10 pair merges; E-151 sampled capped BFS. New here: the full depth-8 edge tally with merge edges separated, giving 0 failing key-kept edges at n = 8 c2. Not in RETRACTIONS (grepped T10-related entries only by identifier; none concern this walk).

## Code changed

New `workshop/rounds/047/experimentalist_deepreplay.py`, a patched copy of `workshop/rounds/019/toolsmith_walk.py` (adds edge/merge/taint tally and a failing-edge list; library untouched). No tests run because no library file changed. Resume was checked only by the count agreement with E-096 (24 316 / 63 221 / 2), as in E-096.

## Next

- Overnight proposal (for OVERNIGHT.md): `.venv/bin/python -u workshop/rounds/047/experimentalist_deepreplay.py 8 --class 2 --depth 9 --budget-sec 540 --ckpt /data/c2.ckpt`, re-run until exit 0 or 2 with "depth 9" printed (about 4-5 slices of 9 minutes; checkpoint grows to above 100 MB, keep outside the repository). Adds one more level where E-147/E-151 failures lived at n = 7 (parent depth 7-8, so depth 9 children are the first with parents at depth 8 for n = 8 c2 already covered here; the point is whether any show up at all).
- The n = 7 classes 1, 2 (where E-151 found failures) with this same tally at the E-151 depth would show whether the failing key-keepers sit on merge paths; that is cheap (minutes) and is the right control for this null: toolsmith or experimentalist next round.
- Skeptic: re-run slice 1 and compare the one-failure state at expansion 14 825 (1 rejection) with the table above.
