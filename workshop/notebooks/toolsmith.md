# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min.
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'until ...; do sleep 5; done'`, <= 100 s per call); `pkill -f` with a pattern in own command line kills own shell; `rm -f $VAR/*` blocked; > 3 jobs on 4 cores inflates timings (seen again 065); a call that times out/exits 144 leaves the background jobs running, check `ps`.

## Round 019-043
- Walk/replay/longsquare/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family.
- `043/toolsmith_revcontrol.py`: reverse control 12/12; 10.3% reverse edges lost (E-149).

## Round 045-053 (T10)
- 047 `merges.py --witness`; 049 meet-in-the-middle; 050/051 all 25 E-151 children joined under J = 0 (E-160, E-163). `canonicalKey` None for bundles > 720 relabelings; cap 5040 changes no verdict.

## Rounds 055-063
- S-1 sizing n = 13-15 (`allRelationLengths`, full n = 15 ~33 min). End(T) iso at every step, J = 0 or not (E-174): only Hom(T,T[-1]) tests the premise. 061 Hom test n = 10 35/35; 063 dim assertion 370 + 35 edges, drop-one control 1 228/1 228.
- Open housekeeping: DEFAULT_CAP 5040 + docstring rewords; relabelling-aware `meetingPoints` (E-164).

## Round 065 (T10 item 1: sizing a J != 0 path at n = 8, 9)
- `065/toolsmith_jplan.py` (wraps `047/experimentalist_deepreplay.py`, parses levels/edges/J != 0). Key-kept edges: n = 7 25 J != 0 of 174 804 (1.4e-4, first at parent depth 7, 110-300 s); n = 8 classes 0-4: 0 of 269 965 (depth 6-9); n = 9 classes 0, 1, 2, 4: 0 of 52 190 (depth <= 7). Expected ~38 at n = 8 at the n = 7 rate.
- Cost: x2.4 expansions per depth; 0.02 / 0.04 / 0.1 s per expansion at n = 7 / 8 / 9. n = 8 c2 depth 9 ~ 20 min, depth 10 ~ 50 min, plus a 6-7 deep return ball: multi-hour sharded, n = 9 hopeless. Not run; no hit.
- Believe: the premise question is not reachable by a bigger ball at n >= 8; the cheap clean test is at n = 7 (is the LNA X reached from a failing child certified equivalent to the seed A?).
- Next: if asked, n = 8 c2 depth 9-10 as OVERNIGHT; the n = 7 X ~ A certification (theorist/skeptic); then the open housekeeping above.
