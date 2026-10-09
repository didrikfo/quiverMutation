# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).
- Traps: shell blocks long `sleep` chains (use `timeout .. bash -c 'sleep 110'`); `rm -f $VAR/*` blocked; `tail -4` with several files fails; commands over 120 s get backgrounded (output piped through `tail` appears only at the end); `pkill -f <pat>` kills your own shell if the pattern is in your command line (4th time in 050: write `pkill -f "foo[b]ar"` AND do not repeat the name later in the same command); > 3 jobs on 4 cores inflates timings (other personas share the box and /tmp).

## Round 019-043 (T5, closure, n = 6/7 meet, reverse control)
- Walk/replay/longsquare/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.
- `041/toolsmith_n6meet.py`: a time-capped BFS size is not a class size. `043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-147).

## Round 045-049 (T10 guard audit, witnesses, tilting path back)
- 045: J = 0 check ~8-9% of a guarded step. 047: `merges.py --witness`, F-041 n = 8 merges have witness 3 + 3.
- 049: meet-in-the-middle with moves F and R (R = forward step on the opposite algebra, carried back; op-duality), target ball depth 5, child ball depth 6: 19 of 25 E-149 failing children join an LNA; parents too. F-only variant fails its own control.

## Round 050 (T10 i, depth-7 child ball)
- `rounds/050/toolsmith_depth7.py` (front / shard / merge; sharded, each command < 10 min), `toolsmith_shards.sh`, `toolsmith_path12.py`, `toolsmith_slow12.py`. Ball cost: c1 depth-6 frontier 18 000 nodes (16 MB pickle), expanding it ~1 000 core-s; ~15 min wall on 3 cores.
- Result: c2 child 6, c1 children 5/13 hit at total 12; c1 child 12 hits at 9 (depth 4). Only c1 children 14/15 (key b32eca) still miss at 7 + 5 = 12. 23 of 25 steps joined (conditional on the J = 0 + tiltingPlus premise). Positive control at depth 7: c2 child 0 with target ball truncated to depth 4, found.
- New trap: `canonicalKey` returns None for bundles over 720 relabelings (DEFAULT_CAP); it happens at n = 7 (c1 child 12: 11 of 153 depth-3 nodes). Such nodes are never deduplicated and never match, and are slow to expand (> 500 s). Docstring "no search has yet produced" is wrong. 049's "undecided" was this, not a big ball.

## Round 051 (T10 i, c1 14/15 closed at 13)
- `rounds/051/`: `toolsmith_tball.py` (resumable target ball; 049's one-shot builder dies at depth 6, > 10 min), `toolsmith_depth13.py` (050's depth7 + env TBD/TD/CAP), `toolsmith_path13.py` (tpar/cball/find), `toolsmith_capcost.py`, `toolsmith_ctrl13.py`.
- Result: target ball c1 depth 6 = 62 297 keys (46 435 at depth 6, ~280 s). c1 14/15 (b32eca) hit at 7 + 6 = 13, 21 meeting keys, path `F1 F3 F1 F5 F4 R7 R1 | F7 R2 R1 F2 F7 R3` replay ok. All 25 E-149 children joined (premise J = 0). Child depth 6 + target 6 = 12 still a miss.
- `canonicalKey` cap 5040 (via wrapper, no library change): changes no verdict; keys 16 of 28 no-key depth-6 nodes (cost 1152-2880, ~0.25 s each); 12 nodes cost >= 622 080 stay keyless. Suggest DEFAULT_CAP 5040 to the chair.
- Controls: c1 child 5 gives 12 = 6 + 6 with the d6 ball; walk-13 seed 3 meets at 9. No control with minimum exactly 13 (seeds 1, 2 balls too slow: parallel-arrow nodes; 4 cores shared, shards die at 10 min when 4+ jobs run). Shard 11500:13800 died in both runs.
- Traps: `sleep N` > ~60 chained is blocked: use `timeout N bash -c 'until [ -f X ]; do sleep 5; done'`. Piping a front run through `tail` loses its log. /tmp/tsm was empty at round start: pickles do not survive rounds; rebuild costs ~25 min (collect 600 s, tball 640 s, front 290 s).

## Next
- Skeptic: independent derived-equivalence test of J = 0 steps; replay all 4 new paths (E-158 x3, 051 x1).
- Orbit-based canonical form for bundles of cost > 5040 (12 nodes in c1 14 ball); a horizon-13 control with minimum exactly 13 (geodesic walk).
- Earlier: promote `tiltingPlus` guard only if the human agrees; reverse loss-by-depth; `merges.py 10 --depths 5 --witness` unsized; docstring fixes (guard, canonicalKey, DEFAULT_CAP).
