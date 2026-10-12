# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots`.
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'while pgrep ..; do sleep 5; done'`, <= 115 s per call; Monitor events fire early and are not reliable for waiting); `rm -f $VAR/*` blocked; > 3 jobs on 4 cores inflates timings; `pkill -f <pat>` kills own shell.

## Round 019-043
- Walk/replay/longsquare/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.
- `043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-149).

## Round 045-053 (T10)
- 047 `merges.py --witness`; 049 meet-in-the-middle; 050/051 all 25 E-151 children joined under J = 0 (E-160, E-163).
- `canonicalKey` None for bundles > 720 relabelings; cap 5040 changes no verdict (suggest DEFAULT_CAP 5040 + docstring rewords + E-160 wrapper test, still open).
- 053 End(T) quiver-level iso on 13/13 E-163 edges.

## Round 055 (S-1 sizing)
- `allRelationLengths(n)` n = 13/14/15: 208 012 / 742 900 / 2 674 440 rows; full n = 15 scan ~33 min on 4 cores; S-1 NOT closed.

## Rounds 057-059 (T10 item 1)
- End(T) iso to the mutation algebra at every step, J = 0 or not (E-174), so only Hom(T,T[-1]) tests the premise. `057/toolsmith_endt2.py` (`symcheck2`, `compare2`, dimcheck from 059: drop-one-relation 125/125 rejected).
- Open: equal-dim wrong algebra at step 12 / c2; DEFAULT_CAP 5040 + docstring rewords; relabelling-aware `meetingPoints` (E-164).

## Round 061 (T10 item 1: Hom(T,T[-1]) at n = 10)
- `061/toolsmith_homn10.py` (`plan | n10 | ctrl PK`): independent `skeptic_tilt.stepTest` on the 35 edges of the 5 E-169 paths: Hom(T,T[-1]) = Hom(T,T[1]) = 0 and Cartan(End T) = child on 35/35 (3 s total). Control: n = 7 c2 failing steps 9/9 Hom(T,T[-1]) = 1, rejected, Cartan still equal. Pickle rebuild c2 = 377 s.
- Believe: group-A merge has no J != 0 edge at the Hom level either; expected from the gate, so this is a cross-check. Power at n = 10 checked (`061/toolsmith_homn10_power.py`): stepTest rejects 312/312 wrong-vertex/direction steps, agrees with tiltingPlus on 700/700; no genuine J != 0 step known at n = 10.
- Next: Hom replay on the other E-163 paths; the open housekeeping above.
