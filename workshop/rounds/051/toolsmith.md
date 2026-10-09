# At n = 7, class 1, under the J = 0 premise, a depth-6 target ball joins c1 children 14 and 15 (key b32eca) to an LNA/dual by a replayed `tiltingPlus` path of total length 13 (7 + 6), so all 25 E-149 failing children are joined; a canonicalKey cap of 5040 changes no verdict

author: toolsmith · round: 051 · kind: result
thread: T10 (i) · bears on: H-015, E-145, E-149, E-155, E-158
scope: n = 7, key class 1 only (c1 children 14 and 15 share the one key b32eca); tilting-only moves F and R (J = 0, `tiltingPlus`, gate, legal, vertex-preserving, class key kept); meet on equal `canonicalKey` between a child ball of depth 7 (one frontier shard of 11500:13800 never finished, not needed for a hit) and the target ball (LNAs + duals) of depth 6 = 62 297 keys, no key missing. Conditional, as in E-155/E-158, on the premise that J = 0 + `tiltingPlus` steps are derived equivalences (E-159 tests 324 edges at Cartan level, generation assumed). A hit is a path of this kind of length 13; a miss would have been a bound only. n = 8 and E-145's own steps not covered.

## Claim

E-158 left the key b32eca (c1 children 14 and 15, parent depth 8, v = 7) open at total length 12. With the target ball extended from depth 5 to depth 6 (62 297 keys: levels 12, 88, 320, 954, 3086, 11402, 46435) and the child ball at depth 7, **there is a path of total length 13: child side `F1 F3 F1 F5 F4 R7 R1` (7 moves), LNA/dual #9 side `F7 R2 R1 F2 F7 R3` (6 moves)**, replayed with fresh move generation, keys equal ("replay ok"). 21 distinct meeting keys at level 7, all total 13; none at child depth 6 (6 + 6 = 12 is a miss, consistent with E-158's 12). Child 15 has the same key, so it is joined by key equality (no separate path printed). With E-158's 23, **all 25 failing children of E-149 are now joined to an LNA/dual by J = 0 `tiltingPlus` paths of length <= 13**, under the premise. Raising the `canonicalKey` cap from 720 to 5040 does not change the outcome (same 21 hits, same levels up to 17 nodes) and is cheap; it keys 16 of 28 no-key nodes of the level-6 frontier, and 12 nodes with cost >= 622 080 stay keyless at any reasonable cap.

It does not say the key guard is right, nor that 13 is the shortest path (13 is the first total at which a hit appears with these two depths; 12 was missed at 6 + 6, 7 + 5 and, by E-158, nothing shorter). It is refuted as a path by any failure of the replay; the premise is the open part.

## Evidence

Sizing first: `toolsmith_depth13.py` levels from E-158 (frontier 18 002 at depth 6), `toolsmith_capcost.py` (cost distribution). Target ball: 62 297 keys, built in 3 resumable stages (round 049's one-shot builder exceeds 10 min at depth 6: it died without a cache, so `toolsmith_tball.py` checkpoints). 0 target nodes without key (level sizes = key counts).

| step | c1 child 14 (b32eca), `canonicalKey` cap 720 | same, cap 5040 |
|---|---|---|
| child levels 0..6 | 1, 10, 56, 229, 982, 4246, 18002 | 1, 10, 56, 229, 982, 4246, 18001 |
| ball nodes through depth 6 | 23 496 | 23 513 |
| best at depth <= 6 (target 6) | none (6 + 6 = 12 miss) | none |
| depth-7 frontier shards expanded | 7 of 8 finished (15 702 + 1 902 nodes; 11500:13800 killed by timeout) | same 7 |
| new depth-7 keys | 69 088 + shard 16100 | 69 151 + shard 16100 |
| no-key children at depth 7 | 239 + 1 | 154 + 0 |
| hit keys, total | 21, all 13 | 21, all 13 |

No-key frontier nodes at depth 6 (cap 720): 28 of 18 002, costs 1152 x12, 1296 x2, 2880 x2 (cap 5040 keys these 16) and 622 080, 21 772 800, 3.3e8, 6.6e8, 7.0e8, 1.9e10, 3.8e13 (12 nodes; no cap helps). `canonicalKey` at cap 5040 on cost 1152-2880 nodes: ~0.25 s each, the same at 40 320. So the cap 5040 costs almost nothing and is the right default for n = 7, but it fixes only the cheap part of the blind spot: the 12 high-cost nodes need a canonical form by orbit, not by enumeration.

Positive controls (same code, same depth):
1. c1 child 5 (key f7abe9, E-158: total 12 = 7 + 5): with the depth-6 target ball the front finds best (12, 6, 6), total 12 at child depth 6, no hit earlier. So the new target ball and the meet code reproduce the known total, split 6 + 6 instead of 7 + 5.
2. A random forward J = 0 walk of length 13 from a class-1 LNA/dual (`toolsmith_ctrl13.py`, seed 3, cost-1 end node): the front finds a meet at (9, 3, 6), i.e. total 9 <= 13 through a depth-6 target key. This checks the depth-6 ball; it is weak (the walk is not geodesic).
3. NOT achieved: a control whose minimum is exactly 13 at the level-7 shard code. Seeds 1 and 2 (walks of length 13) have balls whose depth-6 level did not finish in 10 minutes on a shared box (seed 1: a slow parallel-arrow node at depth 4; seed 2: level 6 beyond 10 min), so they were abandoned unrun. The level-7 shard code is the one that gave E-158's hits at 12 (c2 child 0 control, c1 child 5), and here finds its own hit with a replay-ok path; for a hit, the replay is the control, since the path is checked move by move.

Because the result is a hit, the usual worry about misses (caps, exclusions, shard 11500:13800, the unexpanded node 12415 of E-158) does not apply: they can only hide further hits.

## Reproduction

From the repository root, scratch in `/tmp/tsm` (pickles, not committed); times on a 4-core box, 3 jobs at a time.
```
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100   # run twice (checkpoint, 520 s + 88 s)
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/051/toolsmith_tball.py 1 6     # run twice: 520 s, then 119 s -> tball_c1_d6.pkl, 62 297 keys
timeout 10m .venv/bin/python -u workshop/rounds/051/toolsmith_depth13.py front /tmp/tsm/c1.pkl 1 14 7 c1_14     # 288 s, frontier 18 002 (CAP=5040 for the cap run: 304 s)
NODELIM=100 JOBS=3 workshop/rounds/051/toolsmith_shards.sh c1_14 2300 18002    # 8 shards, 275-600 s each; then: .venv/bin/python workshop/rounds/051/toolsmith_depth13.py merge c1_14
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/051/toolsmith_path13.py /tmp/tsm/c1.pkl 1 14 tpar   # twice (520 s, then ~10 s)
timeout 10m .venv/bin/python -u workshop/rounds/051/toolsmith_path13.py /tmp/tsm/c1.pkl 1 14 cball             # 301 s
SL=0:2300 NODELIM=100 timeout 10m .venv/bin/python -u workshop/rounds/051/toolsmith_path13.py /tmp/tsm/c1.pkl 1 14 find   # 40 s: prints the path and "replay ok total 13"
.venv/bin/python workshop/rounds/051/toolsmith_capcost.py /tmp/tsm/fr_c1_14.pkl      # seconds
```
Control 1: `front /tmp/tsm/c1.pkl 1 5 7 c1_5` (313 s). Control 2: `toolsmith_ctrl13.py /tmp/tsm/ctrl3.pkl 13 3`, then `front /tmp/tsm/ctrl3.pkl 1 0 7 ctrl13_3` (stopped after level 5; the meet appeared at level 3). Logs: `workshop/rounds/051/toolsmith_logs.txt` (5 KB).

## Prior record

E-158 (c1 14/15 open at 12; no-key diagnosis), E-155 (19 of 25), E-153 (no-key nodes at n = 8), E-159 (independent Hom test of the printed paths, Cartan level). The reported target ball of depth 6 and the cap numbers are new; grepped `research/` for "b32eca" and "cap 5040": only E-158 for the key. The "cap fix costs recall only" reading follows `fingerprint.canonicalKey`'s docstring (a None is treated as unseen).

## Code changed

New scripts only, in `workshop/rounds/051/`: `toolsmith_depth13.py` (copy of 050's `toolsmith_depth7.py` with env `TBD` target-file depth, `TD` truncation default 6, `CAP` key cap by wrapping `fingerprint.canonicalKey`), `toolsmith_tball.py` (resumable target ball), `toolsmith_path13.py` (tpar / cball / find path recovery at 13), `toolsmith_ctrl13.py`, `toolsmith_capcost.py`, `toolsmith_shards.sh`, `toolsmith_logs.txt`. No library file touched, so no tests run. The `canonicalKey` docstring still says "no search has yet produced" None; I did not edit it.

## Next

- skeptic: replay the path above and the 3 E-158 paths with the independent Hom(T,T[m]) test (E-159 method); it is the only premise left for T10 (i).
- chair/human: consider `DEFAULT_CAP = 5040` for n = 7 searches (0.25 s per keyed node of cost 1152-2880, no verdict changed); the 12 keyless nodes with cost >= 622 080 need an orbit-based key.
- toolsmith: a horizon-13 control with minimum exactly 13 (needs a geodesic walk; or build the control from the found path: truncate the target ball to depth 5 and search child depth 8, an OVERNIGHT-size job, not proposed).
