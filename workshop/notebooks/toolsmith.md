# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).
- Traps: shell blocks long `sleep` chains (use `timeout .. bash -c 'until ...'`, <= ~110 s per call); `rm -f $VAR/*` blocked; `tail -4` with several files fails (use `-n`); commands over 120 s get backgrounded; `pkill -f <name>` kills your own shell if the name is in the command line; > 3 jobs on 4 cores inflates timings 2-3x.

## Round 019-033 (T5)
- `rounds/019/toolsmith_walk.py`, `022/toolsmith_replay.py`; `QM_CHECK_CARTAN=1` opt-in. `toolsmith_longsquare.py` + tests/test_longsquare.py. `toolsmith_snfresolve.py`.
- Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.

## Round 035-043 (closure, n = 6/7 meet, reverse control)
- `rounds/041/toolsmith_n6meet.py`: 16 hits share 0 keys with the tilting-only LNA side; a time-capped BFS size is not a class size.
- `rounds/043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-147).

## Round 045 (T10 guard audit a)
- `rounds/045/toolsmith_guardaudit.py`: J = 0 check ~8-9% of a guarded step; key-keeping J != 0 steps exist at n = 7 c1/c2 only. Docstring claim too strong.

## Round 047 (T10 iii, witnesses)
- `merges.py --witness` (opt-in, key `witnesses` in JSONL: shortest path + start 0 member / 1 dual); `tests/test_merges_witness.py`. Note `searchFrom` now returns 7 items.
- `rounds/047/toolsmith_witness.py [--all]`: replays `search.meetingPoints` paths (negative steps = opposite algebra) with gate/tiltingPlus/key per edge. Result: the 10 F-041 n = 8 merges have shortest witness 3 + 3 = 6 (not <= 4 in total), all 60 edges J = 0 and key-keeping; B halves replayed forward, inverse edges untested.
- Traps: exec of `scholar_h015.py` clobbers names `ap` (= arrowPaths); relLengths has n - 2 entries (an LNA tuple of n entries fails); F-041 "merges" are `meetingPoints`, not `merges.py` links, so `merges.py` witnesses are untested on real links (no link at n = 9 depth 3).

## Next
- Promote `tiltingPlus` + `tiltingGuard=False` only if the human agrees. Reverse loss-by-depth; the 3 h reverse job is not in OVERNIGHT.md. Deep merges (`merges.py 10 --depths 5 6 --witness`) not sized. 'nokey' counter question from 033 open.
