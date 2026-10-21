# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- Traps: long `sleep` chains blocked (<= 115 s per call, `timeout N bash -c 'until ...'`); `pkill -f` with a pattern in own command line kills own shell; > 3 jobs on 4 cores inflates timings; a timed-out call leaves jobs running, check `ps`.

## Round 019-065 (short)
- Walk/replay/tiltpath/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS/EXPERIMENTS before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control. `canonicalKey` None for bundles > 720 relabelings; cap 5040 changes no verdict. End(T) iso at every step (E-174): only Hom(T,T[-1]) tests the premise.
- Open housekeeping: DEFAULT_CAP 5040 + docstring rewords; relabelling-aware `meetingPoints` (E-164).

## Round 067-075 (T10/T5)
- E-168/E-198 scanned paths for cycles; 70 key-kept failing edges (34 + 36), not capped 25 (E-204). Check whether an EXPERIMENTS count was capped.
- 073: n = 7 class 2, 8 of 9 E-151 failing steps parent ~ child by Hom-verified chains (needs `049/toolsmith_collect.py`, 311 s). `canonicalKey` is LABELLED; no NO-test has power at n = 7. n = 9 F-010 5 + 5 is an OVERNIGHT proposal.
- 075: iso-meet J = 0 joins of the 4 unjoined E-213 classes (child 2,2,2,3; parent 6,6,5,6). Always meet up to sign-aware iso; labelled misses are not evidence.

## Round 078 (chains recorded and Hom-verified; E-211 refusal)
- Believe: all 8 E-216 chains (32 edges) pass the E-212 Hom+Cartan test and iso meets. `078/toolsmith_chains.py R 4 5` (R=0..3, ~85 s) records the move lists (parent pointers) and verifies; it execs the 073 certificate prefix (up to `recs = [r for r`), 075 join prefix, 075 isojoin helpers. Still missing: parent ~ child seed-to-seed link (the chains end at seeds; not printed which).
- Added `CyclicDescendantsError` + `allowCyclic` flag in `procedure.mutateAtVertex` (D(i) cyclic -> refuse), test in `tests/test_cyclic_warning.py`; 9 test files pass (not slow). Not run: slow tests, other files; direct callers on cyclic D(i) now raise.
- Traps: a local function named `search` shadows the `search` module the exec'd preamble uses; `run_side` is defined after the pkl-loading line in the certificate, copy it rather than exec past it.
- Next: print seed indices and verify the seed-to-seed link (seeds of the same LNA class, E-157 style); iso-meet default in `meetingPoints`; record 6 of E-157 with iso meets; F-010 5 + 5 overnight.
