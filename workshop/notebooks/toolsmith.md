# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).
- Traps: shell blocks long `sleep` chains (use `timeout .. bash -c 'until grep -q ...'`, or sleep <= ~110 s); `rm -f $VAR/*` blocked; `tail -4` with several files fails (use `-n`); commands over 120 s get backgrounded, so poll; after a backgrounded cd, use `cd` in each command (cwd was lost once).

## Round 019-033 (T5)
- `rounds/019/toolsmith_walk.py`, `022/toolsmith_replay.py`; `QM_CHECK_CARTAN=1` opt-in. `toolsmith_longsquare.py` + tests/test_longsquare.py.
- `toolsmith_snfresolve.py`: places unresolved LNAs via F-047 profile. Lesson: grep FINDINGS for a derived invariant before sizing.
- 031: doubled-arrow parents gate-admitted, Cartan fails; hand-built parallel algebras need `procedure.toPathAlgebra`.
- 033 `toolsmith_baserate.py`: key test vacuous for E-121; "no hit" needs a positive control in the same family.

## Round 035-038 (closure, n = 6/7)
- `rounds/035/toolsmith_closure.py`: n = 7 classes not closable in 10 min. `rounds/038/toolsmith_n6close.py` (`--plan S --only IDX`, `--budget-hours`), `toolsmith_single.py` (E-133 crash: label dicts keyed by rows of length n-2).
- n = 6 has 4 LNA key classes (2, 24, 26, 32 LNAs); none closes by BFS in minutes (~130-190 seen/s).
- Lesson: meet-in-the-middle with canonical keys; bounded "no hit" must be re-asked at 100x the bound.

## Round 041 (tilting-only meet, E-137 follow-up)
- `rounds/041/toolsmith_n6meet.py`: `--tilting-only --control --revcontrol --reverse --hits --plan --budget-hours`; prints CLOSED / CAP HIT per BFS.
- Positive control (node z with two tilting parents; BFS from each must share z) passes. Reverse search = forward tilting step on opposite algebras; only ~87% of LNA-side tilting edges are inverted that way, so reverse misses are weak.
- All 16 hits: 0 shared with LNA side, forward and reverse. Nothing closes; hit 0 forward passes 5 356 and keeps growing (E-137's 323-501 were a 12 s cap). J != 0 steps are ~0.1% (LNA) / 6% (hit) of gate-admitted steps, so tilting-only barely prunes.
- Lesson: a time-capped BFS size is not a class size; always print the frontier left.

## Next
- Overnight (proposed in the 041 submission): 3 h tilting-only meet with --reverse on hits 0, 4, 13.
- Skeptic: why 13% of tilting edges are not inverted by opposite-step.
- Better than BFS-vs-BFS: an invariant that separates the hit from the class.
- Open from 033: does `_coxeterKeyOrNone` return None on any child (silent edge drops)? The 'nokey' counter exists in the new bfs but is not printed.
