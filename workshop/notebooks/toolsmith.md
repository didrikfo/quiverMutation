# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots`.
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'until grep -q X f; do sleep 5; done'`); `rm -f $VAR/*` blocked; > 3 jobs on 4 cores inflates timings; `pkill -f <pat>` can kill own shell.

## Round 019-043
- Walk/replay/longsquare/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.
- `043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-147).

## Round 045-053 (T10: tilting path back; quiver-level End(T))
- 047 `merges.py --witness`; 049 meet-in-the-middle (moves F, R); 050/051 all 25 E-149 children joined under J = 0 (E-158, E-161).
- `canonicalKey` returns None for bundles > 720 relabelings; cap 5040 changes no verdict (suggest DEFAULT_CAP 5040 in a toolsmith round with docstring rewords + E-158 wrapper test).
- 053 `toolsmith_endt.py`/`_run.py`: End(T) ~ next algebra on 13/13 E-161 edges, but also at 8 decided failing J != 0 steps, so it cannot test the premise. Gaps: parallel-arrow steps, class 2, E-158/E-155 paths.
- Traps: exec-ing slices of other scripts misses module constants; pickles in /tmp/tsm do not survive rounds (`toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100` rebuilds in 513 s).

## Round 055 (S-1 sizing)
- `allRelationLengths(n)` returns a tuple (fully built): n = 13/14/15 = 208 012 / 742 900 / 2 674 440 rows, 1.3/6.2/23 s. `lnaCoxeterKey` ~3 ms/row at n = 15 -> full key scan ~8000 CPU-s, ~33 min on 4 cores in 24 shards.
- Cheaper test exists: one-row witness. Lone 3 (5,6)/(6,5) at n = 15 share a key, head/tail `removeVertex` images differ in key (also n = 17 (6,7)). `055/toolsmith_s1witness.py`, `toolsmith_s1plan.py`. Scan would only add counts; proposed S-1 closure / OVERNIGHT park.
- Lesson: ask whether the predicted failure needs the class scan or only a member of the class.

## Next
- Left approximation test P_v -> add(T/P_v) on the 13 E-161 edges (AI 2.31); orbit canonical form for the 12 keyless nodes; DEFAULT_CAP 5040 + docstrings; relabelling-aware `meetingPoints` (E-162).
