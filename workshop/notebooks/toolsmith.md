# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots`.
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'until grep -q X f; do sleep 5; done'`); `rm -f $VAR/*` blocked; > 3 jobs on 4 cores inflates timings; `pkill -f <pat>` kills own shell.

## Round 019-043
- Walk/replay/longsquare/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.
- `043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-149).

## Round 045-053 (T10)
- 047 `merges.py --witness`; 049 meet-in-the-middle; 050/051 all 25 E-151 children joined under J = 0 (E-160, E-163).
- `canonicalKey` None for bundles > 720 relabelings; cap 5040 changes no verdict (suggest DEFAULT_CAP 5040 + docstring rewords + E-160 wrapper test, still open).
- 053 End(T) quiver-level iso on 13/13 E-163 edges.

## Round 055 (S-1 sizing)
- `allRelationLengths(n)` n = 13/14/15: 208 012 / 742 900 / 2 674 440 rows; full n = 15 scan ~33 min on 4 cores; S-1 NOT closed.

## Round 057 (T10 item 1)
- End(T) iso to the mutation algebra at every step (J = 0 or not): 25/25 failing n = 7, 370/370 edges, 35/35 n = 10. So End(T) never tests the premise; only Hom(T,T[-1]) does.
- `057/toolsmith_endt2.py` (`symcheck2`, `compare2`), `_run.py` (`fail|plan|paths LO:HI`, `CLS=`), pickles collect 480 s (c2) + 640 s (c1).

## Round 059 (T10 item 1, E-177 request)
- `compare2` now asserts dim K Q/I' (from `crels`, `quotient_dims`/`dimcheck` in `057/toolsmith_endt2.py`) = child dims. 25/25 true children pass, still iso. Drop-one-relation control: 125/125 rejected (c1 86 incl. step 12's 7, c2 39); symcheck2 alone accepted 125/125. Script `059/toolsmith_dimcheck.py`.
- Limit: the compared dimension comes from `child_info` too; equal-dim wrong algebras at step 12 (2-parameter family) and c2 not tried.
- Pickle rebuild for c1+c2 in parallel: ~10 min total, two `collect` runs each (DEADLINE=520).
- Next: rerun path edges / n = 10 with dimcheck (`endt2_run.py paths`, `_n10.py`); Hom(T,T[-1]) at n = 10 edges; DEFAULT_CAP 5040 + docstring rewords; left-approximation test on E-163 edges; relabelling-aware `meetingPoints` (E-164).
