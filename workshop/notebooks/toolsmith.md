# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).
- Traps: shell blocks long `sleep` chains (use `timeout .. bash -c 'until grep -q ...'`, sleep <= ~110 s); `rm -f $VAR/*` blocked; `tail -4` with several files fails (use `-n`); commands over 120 s get backgrounded, so poll; after a backgrounded cd, use `cd` in each command.

## Round 019-033 (T5)
- `rounds/019/toolsmith_walk.py`, `022/toolsmith_replay.py`; `QM_CHECK_CARTAN=1` opt-in. `toolsmith_longsquare.py` + tests/test_longsquare.py. `toolsmith_snfresolve.py` (F-047 profile).
- Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family (033 baserate).

## Round 035-041 (closure, n = 6/7 meet)
- `rounds/038/toolsmith_n6close.py` (`--plan S --only IDX`), `toolsmith_single.py`. n = 6 has 4 LNA key classes (2, 24, 26, 32 LNAs); none closes by BFS in minutes (~130-190 seen/s).
- `rounds/041/toolsmith_n6meet.py`: `--tilting-only --control --revcontrol --reverse --hits --plan --budget-hours`; prints CLOSED / CAP HIT. All 16 hits share 0 keys with the tilting-only LNA side, fwd and rev; nothing closes. J != 0 steps ~0.1% (LNA) / 6% (hit).
- Lesson: a time-capped BFS size is not a class size; print the frontier left. Meet-in-the-middle with canonical keys.

## Round 043 (reverse control, E-142 follow-up)
- `rounds/043/toolsmith_revcontrol.py` (`--diag`, `--control --depth k`, `--budget-hours`). Reverse control passes: forward path A->B->C of length 2, 3, 4 from a start LNA is found backward from C at exactly depth k, 12/12 each (reverse sets 10-120 nodes, closed).
- Lost reverse edges (10.3% of 1500, first-parent edges): never a filter. Opposite(B) mutated at v gives A2 != A (same class key, same arrow count, 72/154 already in LNA BFS); no vertex gives A. Reverse graph is not the transpose of the forward graph but stays inside the class.
- Caveat: controls have a forward witness; short paths survive 24 edges at ~8% odds, so loss may not be uniform in depth (not tallied).

## Next
- Chair adds the 3 h job (`--tilting-only --reverse --hits 0,4,13`) to OVERNIGHT.md; now unblocked.
- Skeptic/theorist: why mutation of opposite(B) at v is not the inverse in ~10% (dump the 154 edges; shape lead: 6-7 arrows, commutative square into B).
- Loss by depth tally; a reverse-only control.
- Better than BFS-vs-BFS: an invariant separating the hit from the class.
- Open from 033: does `_coxeterKeyOrNone` drop children silently? ('nokey' counter exists, still not printed in the 041 script; the 043 script keeps it in stats but does not print it).
