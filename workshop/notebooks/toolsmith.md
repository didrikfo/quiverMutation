# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min.
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'until ...; do sleep 5; done'`, <= 100 s per call); `pkill -f` with a pattern in own command line kills own shell; `rm -f $VAR/*` blocked; > 3 jobs on 4 cores inflates timings; a call that times out/exits 144 leaves jobs running, check `ps`.

## Round 019-043
- Walk/replay/longsquare/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.
- `043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-149).

## Round 045-053 (T10)
- 047 `merges.py --witness`; 049 meet-in-the-middle; 050/051 all 25 E-151 children joined under J = 0 (E-160, E-163). `canonicalKey` None for bundles > 720 relabelings; cap 5040 changes no verdict.

## Rounds 055-065
- S-1 sizing n = 13-15. End(T) iso at every step, J = 0 or not (E-174): only Hom(T,T[-1]) tests the premise. 065 (E-192): J != 0 join search at n = 8, 9 not reachable (0 hits in 181k edges; multi-hour).
- Open housekeeping: DEFAULT_CAP 5040 + docstring rewords; relabelling-aware `meetingPoints` (E-164).

## Round 067 (T10/T5 item 1: cyclic parents)
- LESSON: grep EXPERIMENTS for the property before scanning. E-168 (054) had already scanned all 44 n = 7 paths (370 edges) for loops AND oriented cycles: 0. STATE/E-193 kept listing it as open. I only added n = 10 E-169 (35 edges, 0) via `067/toolsmith_cycscan.py`.
- Added `CyclicQuiverWarning`/`hasOrientedCycle`/`warnIfCyclic` in `arrowPaths.py`, called in `cartanMatrix`, `arrowPaths.isInIdeal`, `procedure._minimalGenerators`; test `tests/test_cyclic_warning.py`. `relationAlgebra` (vertex model) unguarded. The search does meet cyclic quivers (test_cycles warns).
- Believe: printed paths are in the acyclic setting; the live gap is cyclic children turned parents inside key-kept balls (uncounted).
- Next: count cyclic parents in E-192 balls if asked; DEFAULT_CAP/docstring housekeeping; n = 8 c2 depth 9-10 only as OVERNIGHT.
