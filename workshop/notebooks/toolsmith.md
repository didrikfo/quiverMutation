# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'until grep -q X f; do sleep 5; done'`, < 120 s per call); `rm -f $VAR/*` blocked; commands > 120 s get backgrounded; `pkill -f <pat>` kills own shell if the pattern is in your command line; > 3 jobs on 4 cores inflates timings (other personas share the box).

## Round 019-043 (T5, closure, n = 6/7 meet, reverse control)
- Walk/replay/longsquare/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.
- `043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-147).

## Round 045-051 (T10 guard audit, tilting path back)
- 047: `merges.py --witness`. 049: meet-in-the-middle, moves F and R (R = forward step on the opposite algebra, carried back). 050/051: depth-7 child ball, depth-6 target ball; all 25 E-149 children joined under J = 0 (E-158, E-161).
- `canonicalKey` returns None for bundles over 720 relabelings (n = 7 parallel arrows); cap 5040 changes no verdict, keys 16 of 28 no-key nodes; 12 nodes stay keyless. Suggest DEFAULT_CAP 5040 (chair: not yet, do in a toolsmith round with docstring rewords + E-158 wrapper test).
- Pickles in /tmp/tsm do not survive rounds; `toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100` rebuilds c1 in 513 s (one run, DEADLINE=520).

## Round 053 (T10 item 1: quiver-level End(T))
- `rounds/053/toolsmith_endt.py` (class `TiltEnd`: Hom_K(T_i,T_j) basis, composition, exact over Q) and `toolsmith_endt_run.py` (modes selftest | path13 | tilt | perturb | fail). Radical filtration -> arrows; relation check by Groebner (arrow scalars + radical-square corrections), 2 s per mode.
- Result: End(T) ~ next algebra on 13 of 13 E-161 edges (labelled iso, no parallel arrows). BUT it also agrees at the 8 decided of 16 failing J != 0 steps: End(T) of the 2-term complex is the mutation algebra whether or not T is tilting. So this comparison cannot test the premise; only Hom(T,T[-1]) = 0 (E-159) and generation do. Perturbed relations rejected 22/22 (monomial) -> the decision procedure has power.
- Gaps: parallel-arrow steps undecided (needs matrix-valued arrow identification); no wrong-algebra control from the walk; class 2, E-158 paths, E-155 paths not run; no vertex permutations.
- Trap: exec-ing slices of other scripts (`skeptic_tilt.py`) misses module constants like P; define them.

## Next
- Extend End(T) to class 2, the 3 E-158 paths, 25 E-155 paths (10 min to rebuild pickle each class) -- low value given the negative half; better: theorist's generation / AI 2.31 hypothesis check on the same 13 edges (is P_v -> add(T/P_v) a left approximation?), which I could compute (left add(T/P_v)-approximation test).
- Orbit canonical form for the 12 high-cost keyless nodes; horizon-13 control with minimum exactly 13; relabelling-aware `meetingPoints` (E-162); DEFAULT_CAP + docstring rewords; group-A witness path.
