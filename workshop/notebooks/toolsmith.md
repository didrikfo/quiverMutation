# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots`.
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'until grep -q X f; do sleep 5; done'`); `rm -f $VAR/*` blocked; > 3 jobs on 4 cores inflates timings; `pkill -f <pat>` kills own shell (happened again r057).

## Round 019-043
- Walk/replay/longsquare/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.
- `043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-149).

## Round 045-053 (T10)
- 047 `merges.py --witness`; 049 meet-in-the-middle (moves F, R); 050/051 all 25 E-151 children joined under J = 0 (E-160, E-163).
- `canonicalKey` returns None for bundles > 720 relabelings; cap 5040 changes no verdict (suggest DEFAULT_CAP 5040 + docstring rewords + E-160 wrapper test, still open).
- 053 `toolsmith_endt.py`/`_run.py`: End(T) quiver-level iso on 13/13 E-163 edges, also at 8 decided failing steps.

## Round 055 (S-1 sizing)
- `allRelationLengths(n)` n = 13/14/15: 208 012 / 742 900 / 2 674 440 rows. Full n = 15 key scan ~33 min on 4 cores; one-row witness exists (`055/toolsmith_s1witness.py`). S-1 NOT closed.

## Round 057 (T10 item 1)
- Believe: End(T) of the 2-term complex is iso to the mutation algebra at EVERY step, J = 0 or not (25/25 failing n = 7 steps incl. 8 parallel-arrow; 370/370 edges of E-157/E-160/E-163 paths; 35/35 n = 10 E-169 edges). So End(T) never tests the premise; only Hom(T,T[-1]) does. Lemma candidate for the theorist.
- Did: `057/toolsmith_endt2.py` (`symcheck2`: matrix-valued arrow identification, one Rabinowitsch z per block; z * product of dets hung sympy on a 3x3 + 2x2 block), `_run.py` (`fail|plan|paths LO:HI`, `CLS=` env), `_n10.py`, `_ctrl.py`, `_keys.py`. Whole job < 10 min; collect pickles 480 s (c2) + 640 s (c1).
- Weak: relation-level power is only the monomialisation control (9/9); no equal-dims wrong-algebra control.
- Next: Hom(T,T[-1]) at the n = 10 E-169 edges (skeptic `stepTest`); DEFAULT_CAP 5040 + docstrings; left-approximation test (AI 2.31) on the E-163 edges; relabelling-aware `meetingPoints` (E-164); orbit canonical form for keyless nodes.
