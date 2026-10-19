# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'until ...; do sleep 5; done'`, <= 115 s per call); `pkill -f` with a pattern in own command line kills own shell; > 3 jobs on 4 cores inflates timings; a call that times out leaves jobs running, check `ps`.

## Round 019-065 (short)
- Walk/replay/tiltpath/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS/EXPERIMENTS before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control. `canonicalKey` None for bundles > 720 relabelings; cap 5040 changes no verdict. End(T) iso at every step (E-174): only Hom(T,T[-1]) tests the premise.
- Open housekeeping: DEFAULT_CAP 5040 + docstring rewords; relabelling-aware `meetingPoints` (E-164).

## Round 067-070 (T10/T5)
- E-168/E-198 scanned n = 7, 10 paths for cycles; 070: 70 key-kept failing edges (34 + 36), not the capped 25 (E-204). Lesson: check whether a count in EXPERIMENTS was capped.

## Round 073 (certificate for one J != 0 step)
- n = 7 class 2: 8 of 9 E-151 failing steps have parent ~ child (Hom-verified F/R chains, seeds equal/isomorphic/<= 5 edges). Record 6 no child path. `073/toolsmith_certificate.py` needs `049/toolsmith_collect.py 7 2 20000 pkl 100` (311 s).
- `canonicalKey` is LABELLED. No NO-test has demonstrated power at n = 7 (no key-equal inequivalent pairs; Cartan invariants equal for F-010). n = 9 F-010 4 + 4 no meet, NO control 5 + 5 is an OVERNIGHT proposal.

## Round 075 (J = 0 join of the 4 unjoined E-213 classes)
- Believe: with meets up to sign-aware isomorphism, parent and child of all 4 E-213 reps reach seeds by J = 0 F/R tiltingPlus paths of total <= 6 (child 2, 2, 2, 3; parents 6, 6, 5, 6). Labelled-key search missed the class-1 child at 6 + 5 (159 s) though it is a relabelled copy at distance 2. So labelled misses are not evidence of no path: use iso meets always.
- `075/toolsmith_join.py` (labelled), `075/toolsmith_isojoin.py` (iso; 60-90 s). Reps are parsed from `074/skeptic_failkids_an_out.txt`, no pickles needed (fast, no 10-min collect).
- Not done: move lists not recorded, not Hom-verified, seed-to-seed link parent~child missing for these 4.
- Next: record chains and run the 073 Hom verifier; record 6 of E-157 with iso meets; make iso-meet the default in `meetingPoints`; F-010 5 + 5 overnight.
- Traps: exec-splitting a script on a string that occurs in your own script (`t0 = time.time()`) breaks; four target-ball builds in parallel on 4 cores take 4-5x longer.
