# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).
- Traps: shell blocks `sleep N` chains (use `timeout .. bash -c 'until grep -q ...'` or Monitor); `rm -f $VAR/*` blocked; `tail -4` with several files fails (use `-n`).

## Round 019-026 (T5)
- `rounds/019/toolsmith_walk.py`, `rounds/022/toolsmith_replay.py`; `QM_CHECK_CARTAN=1` opt-in (15 ms/step). `toolsmith_longsquare.py` + tests/test_longsquare.py. n = 9 c0: 14 exp/s, may not close.

## Round 029-033
- `toolsmith_snfresolve.py`: places unresolved LNAs via F-047 profile. Lesson: grep FINDINGS for a derived invariant before sizing.
- 031: doubled-arrow parents are gate-admitted, Cartan fails; hand-built parallel algebras need `procedure.toPathAlgebra`. Preamble hack (exec rounds/023 up to "if a.hand:") overwrites `_argv`.
- 033 `toolsmith_baserate.py`: key test vacuous for E-121 (walk rows carry key by construction; layered family 0/2704). Lesson: "no hit" needs a positive control in the same family.

## Round 035 (n = 7 closure sizing)
- `rounds/035/toolsmith_closure.py plan|run CLASS`: E-124's key-preserving BFS with per-level stats; imports `gen` by exec'ing rounds/033/experimentalist_bothdie.py with argv mode 'none'.
- Class idx 1 (12 LNAs, 38 targets): 48-89 exp/s, frontier ratio 2.2 -> 2.5 (not falling), 29k seen at level 8 (240 s). Class idx 3 (58 LNAs, 6 targets): ~80 exp/s, ratio 2.6 -> 1.9, 40k seen at level 7. 0/44 hits. Not closable in 10 min; class 1 > 10 h by extrapolation.
- Lesson: BFS closure of a key class has no a-priori size; measure the frontier ratio first, sharding by first mutation does not help (shared seen set).

## Next
- Reverse search from the 44 targets toward an LNA (priority by arrows / dims), size with a budget.
- Overnight: `closure.py run 1 --budget-hours 8`, `run 3 --budget-hours 4` (proposed in submission).
- Does `_coxeterKeyOrNone` return None on any child? (could silently drop BFS edges)
- Still open from 033: match 42 non-c0 children to E-114 rejects; positive control family for key test.
