# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).
- Traps: shell blocks long `sleep` chains (use `timeout .. bash -c 'until ...'`, <= ~110 s per call); `rm -f $VAR/*` blocked; `tail -4` with several files fails (use `-n`); commands over 120 s get backgrounded; `pkill -f <name>` kills your own shell if the name is in the command line; with more jobs than cores timings inflate 2-3x (run <= 3 jobs on 4 cores).

## Round 019-033 (T5)
- `rounds/019/toolsmith_walk.py`, `022/toolsmith_replay.py`; `QM_CHECK_CARTAN=1` opt-in. `toolsmith_longsquare.py` + tests/test_longsquare.py. `toolsmith_snfresolve.py`.
- Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.

## Round 035-043 (closure, n = 6/7 meet, reverse control)
- `rounds/041/toolsmith_n6meet.py` (`--tilting-only --control --revcontrol --reverse --hits --plan`): 16 hits share 0 keys with the tilting-only LNA side; no BFS closes. A time-capped BFS size is not a class size.
- `rounds/043/toolsmith_revcontrol.py`: reverse control passes 12/12 (depth 2-4); 10.3% reverse edges lost, not a filter (opposite(B) mutated at v gives another same-key algebra).

## Round 045 (T10 guard audit a)
- `rounds/045/toolsmith_guardaudit.py n cls maxnodes A|J|T`: key-guarded BFS, times gate/step/J/tiltingPlus per call, tallies (J, tp, key kept).
- Believe: J = 0 check costs ~0.4 ms (n=7) / 0.8 ms (n=8) per step, about 8-9% of a guarded step (5 / 9.6 ms); placed right after the gate at search.py line 435, before quiverMutationAtVertex; J == 0 iff tiltingPlus on ~250k steps. Key-keeping J != 0 (tp-failing) steps: 16 (c1), 9 (c2) of ~80k at n=7 after 20k expansions; 0 at n=8 c0/c1 (1200 exp.). Docstring claim is too strong. perI/tiltingPlus are workshop-only; promotion is the human's call (q3).
- Not done: n = 8 c2, deeper n = 8, checking children against a derived-class invariant, a library `tiltingGuard` option.

## Next
- If chair approves: promote `tiltingPlus` to library + `tiltingGuard=False` kwarg + test; fix docstring.
- Loss-by-depth tally for the reverse search; the 3 h reverse job (`--tilting-only --reverse --hits 0,4,13`) is not in OVERNIGHT.md.
- Open from 033: does `_coxeterKeyOrNone` drop children silently ('nokey' counter not printed).
