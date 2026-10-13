# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots`.
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'until [ -f X ]; do sleep 5; done'`, <= 110 s per call, else it backgrounds); `pkill -f`/`pgrep -f | kill` with a pattern in own command line kills own shell (kill by pid from `ps -eo pid,etimes,cmd`); `rm -f $VAR/*` blocked; > 3 jobs on 4 cores inflates timings.

## Round 019-043
- Walk/replay/longsquare/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.
- `043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-149).

## Round 045-053 (T10)
- 047 `merges.py --witness`; 049 meet-in-the-middle; 050/051 all 25 E-151 children joined under J = 0 (E-160, E-163).
- `canonicalKey` None for bundles > 720 relabelings; cap 5040 changes no verdict (suggest DEFAULT_CAP 5040 + docstring rewords + E-160 wrapper test, still open).
- 053 End(T) quiver-level iso on 13/13 E-163 edges.

## Round 055 (S-1 sizing)
- `allRelationLengths(n)` n = 13/14/15: 208 012 / 742 900 / 2 674 440 rows; full n = 15 scan ~33 min on 4 cores; S-1 NOT closed.

## Rounds 057-061 (T10 item 1)
- End(T) iso to the mutation algebra at every step, J = 0 or not (E-174), so only Hom(T,T[-1]) tests the premise. `057/toolsmith_endt2.py` (`symcheck2`, `compare2`, dimcheck from 059).
- 061 `toolsmith_homn10.py`: Hom(T,T[+-1]) = 0 and Cartan equal on 35/35 n = 10 edges; stepTest power checked (312/312 wrong steps rejected).
- Open: DEFAULT_CAP 5040 + docstring rewords; relabelling-aware `meetingPoints` (E-164).

## Round 063 (E-186 open item)
- `063/toolsmith_dimedges.py n7 LO:HI | n10`: dim K Q/I' = dim End(T) = child dims on 214 + 156 = 370 n = 7 path edges and 35 n = 10 edges, 0 failures, none unevaluable. `toolsmith_dimedges_ctrl.py`: 1 228 single-relation drops on those edges all rejected (none redundant).
- Believe: the dimension assertion is a clean but weak check (only catches wrong dimension); nothing new about J. Pickle rebuilds: c1 ~9 min, c2 ~7 min (collect resumes; run both in parallel, not twice).
- Next: housekeeping (cap 5040, docstring rewords, relabelling-aware meetingPoints); edge-level equal-dim wrong-ideal control only if someone asks.
