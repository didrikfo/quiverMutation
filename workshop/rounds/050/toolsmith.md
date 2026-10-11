# At n = 7 (classes 1, 2), under the J = 0 premise, a depth-7 child ball joins c1 child 12 (undecided in E-157) and 3 of the 5 E-157 misses (c2 6, c1 5 and 13) to an LNA/dual by `tiltingPlus` paths of length <= 12 (c1 13 by key equality with 5, no path printed); c1 14 and 15 (one key, b32eca) are open at 12, no-key nodes unmatched

author: toolsmith · round: 050 · kind: result
thread: T10 (i) · bears on: H-015, E-147, E-151, E-154, E-157
scope: n = 7, key classes 1 and 2 (child balls: nodes without a canonical key are expanded, never matched; target balls have none), the 25 E-151 children rebuilt (same collector, 20 000 expansions each, now resumable); tilting-only moves F and R (J = 0, `tiltingPlus`, gate, legal, vertex-preserving, class key kept), meet on equal `canonicalKey` with the depth-5 target ball of the class (LNAs + duals). Conditional, as in E-157, on the premise that J = 0 + `tiltingPlus` steps are derived equivalences (no independent test exists; the skeptic has that item). Total path length 12 is the horizon; n = 8 and E-147's own step list not covered. Nodes without a canonical key (see below) can be expanded but never matched or deduplicated.

## Response to referee

Verdict minor revision; items 1-5 done, item 6 deferred.
1. Count fixed (title, Claim): the 4 newly joined children are 3 of the 5 E-157 misses (c2 6, c1 5, c1 13) plus child 12, which E-157 left undecided and is not a miss. 19 + 4 = 23 of 25 unchanged.
2. No-key counts, from the existing pickles by `toolsmith_nokey_count.py c1_14` (seconds): c1 14 child ball, level size minus keys entered in `seen`: 30 no-key nodes at depths 1-6 in total (levels 1+10+56+229+982+4246+18002 = 23 526 nodes, 23 496 keys; the pickle stores `seen` as a set, so no per-depth split without a 5-min rebuild), plus 343 at depth 7. Target balls (depth <= 5): 0 no-key nodes in both classes (level sizes 12, 88, 320, 954, 3086, 11402 equal the per-depth key counts; c2 likewise 14, 104, 426, 1492, 5172, 18244). So the target side has no missing partner; the loss is one-sided, child side only. The 30 are expanded but not deduplicated, so the child level sizes are node counts, not distinct-state counts, by at most 30 plus their descendants.
3. Table: c1 5 marked "stopped at first hits (6 900 of 18 034)"; c1 13 marked "key equality with 5 only, no path printed".
4. Premise now in title and first line of Claim.
5. E-155 cited (Prior record): n = 8, two depth-8 children without key; the new part is scale (11/153, 100/659, 343).
6. Not done this sitting: a control through the time-limit plus no-key route needs a new ball with a known hit by another route; no such case is built. Child 12 at depth 4 (100 no-key children, hit at 9) stays the nearest. Scope of the claim is unchanged by it: the c1 14 miss is a lower-bound statement with those exclusions.

## Claim

E-157 left 5 children that missed at child depth 6 + target depth 5 = 11 (c2 child 6; c1 children 5, 13, 14, 15; 3 + 1 distinct keys) and c1 child 12 undecided. With the child ball extended to depth 7 (total horizon 12), **under the J = 0 premise, c2 child 6, c1 children 5 and 13 (same key f7abe9; 13 by key equality with 5, no path printed for 13) are each joined to an LNA/dual by a J = 0 `tiltingPlus` path of total length 12** (7 child-side moves + 5 target-side), and **c1 child 12 by a path of total length 9** (4 + 5; E-157 had completed only depth 3 for it). Both new paths replay with fresh move generation and equal keys ("replay ok"). So **23 of the 25 children (19 from E-157 plus these 4) are now joined; 2 children (c1 14 and 15, one canonical key b32eca, v = 7, parent depth 8) still miss: no path of length <= 12 of this kind** (child ball depth 7 = 100 649 nodes, exhaustive except for the exclusions listed under Evidence; target ball depth 5 = 15 862 keys). A miss is a bound, not evidence of being outside the class; refuted as a bound by any path of length <= 12 (e.g. a depth-8 child ball would test 8 + 5 = 13 and cost ~4x more).

It does not say the key-guard is right: it says that, under the premise, the guard's `tiltingPlus` failures at n = 7 c1, c2 are steps between members of one class in 23 of 25 cases, and 2 cases are open at horizon 12.

## Evidence

Sharded run (every command <= 10 min; `toolsmith_depth7.py front` builds the ball to depth 6 and pickles the frontier, `shard lo:hi` expands a slice of the depth-6 frontier and checks each new key against the target ball, `merge` combines). No `--plan` flag exists for this script; sizes were taken from E-157's logs (levels printed at every stage).

| child (key) | child-ball levels 0..6 | new depth-7 keys | hit | total |
|---|---|---|---|---|
| c2 child 6 (0d3498) | 1, 8, 34, 121, 436, 1569, 5703 | 20 667 | yes (all hits total 12) | 12 |
| c1 child 5 / 13 (f7abe9) | 1, 10, 56, 229, 982, 4237, 18034 | 29 368 in 3 of 8 shards (stopped at first hits, 6 900 of 18 034 frontier nodes; ball size a lower bound) | yes; c1 13 joined by key equality with 5 only, no path printed | 12 |
| c1 child 14 / 15 (b32eca) | 1, 10, 56, 229, 982, 4246, 18002 | 77 153 (ball 100 649 nodes) | **no** | miss at 12 |
| c1 child 12 (dc45e9) | 1, 6, 32, 153 (depth 4 done) | 659 | yes (3 hit keys) | 9 |
| positive control: c2 child 0 (1c732e), target ball truncated to depth 4 | 1, 8, 34, 121, 436, 1569, 5703 | 20 665 | yes, 11 keys | 11 |

Positive control (same machinery, same depth 7): c2 child 0 is known (E-157) to meet at 6 + 5 = 11. With the target ball truncated to depth <= 4 (`TD=4`) that path can only be found at child depth 7 + 4; the sharded depth-7 code finds it (11 hit keys, total 11, nothing earlier). So the level-7 code path has demonstrated power at the depth in question; the c1 child 5 run (same class, same pickle, same shard script, same level sizes as the missing c1 14) is a second in-class positive. The missing child 14 differs from 5 only in the pair (key, v = 7, parent depth 8); levels 0..6 are equal up to 0.1%.

Printed paths (child side then LNA/dual side, replayed ok; F v = forward step at v, R v = forward step at v of the opposite algebra, carried back):
- c2 child 6, total 12: `F1 F7 R3 R5 R5 F4 R6` | LNA/dual #0: `F1 F1 F2 R6 R5`.
- c1 child 5 (= 13), total 12: `F1 F3 F5 F4 R7 F3 R1` | LNA/dual #9: `F7 R2 R1 F2 F7`.
- c1 child 12, total 9: `R7 R4 R7 R1` | LNA/dual #6: `F1 F5 F4 R7 F3`.
The tp and key counters are again 0 in every shard: with J = 0 the search graph is gate + legal + J = 0 (the J filter fires, 43-357 times per shard).

Exclusions in the miss (honest list): (a) c1 14: 1 of 18 002 depth-6 nodes (index 12415, a node with parallel arrows) did not finish in 520 s and is unexpanded; 2 others needed 216-231 s on their own and were completed. (b) **Nodes without a canonical key.** `fingerprint.canonicalKey` returns None when a parallel-arrow bundle needs more than 720 relabelings (DEFAULT_CAP). Its docstring says no search has produced this; it occurs here: in the c1 child-12 ball 11 of 153 depth-3 nodes and 100 of 659 depth-4 children have no key, and in the c1 14 ball 343 depth-7 children have none. Such nodes are expanded (not deduplicated) but cannot meet the target ball, so a path whose meeting node has no key is invisible. c1 child 12 was undecidable in 049 because these nodes make the child ball blow up (a node can take > 500 s: 3 of 153 depth-3 nodes needed > 100-150 s), not because it is large: its 4-level ball is 839 nodes. (c) Depth-7 for c1 5/13 was stopped after 6 900 of 18 034 frontier nodes because 9 hit keys had appeared; the ball size for it is therefore a lower bound. The E-157 hits and the 25 parents are unchanged.

Counts after this round: c1 11 + 3 (5, 13, 12) = 14 of 16 steps; c2 8 + 1 = 9 of 9; total 23 of 25 steps (hit keys: 15 + 3 = 18 of 19 distinct; 1 distinct key, b32eca, open).

## Reproduction

All from the repository root, `/tmp/tsm` as scratch (pickles of ~16 MB each, not committed). Times are on a 4-core box shared with another job.
```
timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100   # 370 s (also c2: 274 s); resumable via DEADLINE, did not need it; or toolsmith_run.sh
timeout 10m .venv/bin/python workshop/rounds/049/toolsmith_tiltpath.py /tmp/tsm/none.pkl 1 5 4 400000 ctrl0:9   # target ball c1 depth 5 cache, ~150 s (c2: 196 s)
timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_depth7.py front /tmp/tsm/c1.pkl 1 14 7 c1_14    # 294 s: ball to depth 6, pickles the frontier (18 002 nodes, 16 MB)
workshop/rounds/050/toolsmith_shards.sh c1_14 2300 18002          # 8 shards, 3 at a time, each under timeout 10m: 130-420 s per shard; NODELIM=100 (per-node limit) for the slow ones
.venv/bin/python workshop/rounds/050/toolsmith_depth7.py merge c1_14
TD=4 ... front /tmp/tsm/c2.pkl 2 0 7 c2_ctrl0 ; TD=4 ... shard c2_ctrl0 0:3000 / 3000:6000   # control: 95 s per stage
timeout 10m .venv/bin/python -u workshop/rounds/049/toolsmith_paths.py /tmp/tsm/c2.pkl 2 7 buildball  # 199 s; then SLICE=6:7 ... toolsmith_paths.py ... 2 7 paths (220 s)
DEP=6 NODELIM=60 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_path12.py /tmp/tsm/c1.pkl 1 5   # c1 child 5 path (after buildball for c1, 138 s)
NODELIM=100 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_path12.py /tmp/tsm/c1.pkl 1 12        # child 12 path, 470 s (DEP default 3)
```
Merged and printed logs: `workshop/rounds/050/toolsmith_depth7_logs.txt` (5 KB). `toolsmith_slow12.py` times the child-12 nodes and counts None-key nodes (the diagnosis above).

## Prior record

E-157 (19 hits, 5 misses at 6 + 5, child 12 undecided), E-158, E-151, E-154. `canonicalKey` returning None is noted in E-155/E-131-style entries ("two depth-8 children have no canonical key", "Warning: can return None") and E-155 records the phenomenon at n = 8 (two depth-8 children without key, never entered in `seen`); what is new is the scale at n = 7 (11 of 153, 100 of 659, 343 at depth 7, 30 at depths 1-6 of c1 14) and that it is child-side only. I found no entry saying that a parallel-arrow ball at n = 7 loses matchable nodes; the docstring line "no search has yet produced" is therefore wrong. Grepped `research/` for "no canonical key", "nokey": only those. The "in the class" reading remains conditional on the premise, as in E-157.

## Code changed

New scripts only, in `workshop/rounds/050/`: `toolsmith_collect.py` (049's collector with checkpoint/resume; not exercised, runs finished first), `toolsmith_run.sh`, `toolsmith_depth7.py` (front / shard / merge), `toolsmith_shards.sh`, `toolsmith_path12.py` (path recovery when the ball has no-key nodes), `toolsmith_slow12.py`, `toolsmith_nokey_count.py` (response item 2), `toolsmith_depth7_logs.txt`. No library file touched; no tests run (none apply).

## Next

- skeptic: independent derived-equivalence test of one J = 0 step (still the only premise); replay the three new paths (cheap: no ball build).
- toolsmith (OVERNIGHT proposal): c1 14 / 15 only: 8 + 5 = 13 needs ~300 000 child nodes (about 4x the 100 649 here, ~2 h on 4 cores) or a target ball of depth 6 (57 000 keys, ~20 min once, shardable) with child depth 7 for horizon 13; both bounded by the unexpanded node 12415. A canonicalKey with a larger cap (7! = 5 040 relabelings) would fix the no-key blind spot at a cost per node I did not measure.
- chair: the `canonicalKey` docstring ("no search has yet produced" None) and the T10 (iv) docstring reword wait on the premise test.
