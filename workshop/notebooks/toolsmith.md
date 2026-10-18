# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- Traps: long `sleep` chains blocked (use `timeout N bash -c 'until ...; do sleep 5; done'`, <= 100 s per call); `pkill -f` with a pattern in own command line kills own shell; > 3 jobs on 4 cores inflates timings; a call that times out/exits 144 leaves jobs running, check `ps`.

## Round 019-065 (short)
- Walk/replay/tiltpath/snfresolve scripts; `QM_CHECK_CARTAN=1` opt-in. Lessons: grep FINDINGS/EXPERIMENTS for a derived invariant before sizing; doubled-arrow parents need `procedure.toPathAlgebra`; "no hit" needs a positive control in the same family. `canonicalKey` None for bundles > 720 relabelings; cap 5040 changes no verdict. End(T) iso at every step, J = 0 or not (E-174): only Hom(T,T[-1]) tests the premise.
- Open housekeeping: DEFAULT_CAP 5040 + docstring rewords; relabelling-aware `meetingPoints` (E-164).

## Round 067-070 (T10/T5 item 1)
- E-168 had already scanned n = 7 paths for cycles; E-198 added n = 10 and `CyclicQuiverWarning`. 070: key-guarded n = 7 balls (c1, c2 depth 8) and n = 8 c0/c2 hold no cyclic child; 70 key-kept failing edges (34 + 36), not the capped 25 (E-204).
- Lesson: check whether a count in EXPERIMENTS was capped/checkpointed.

## Round 073 (T10/T5 item 1: certificate for one J != 0 step)
- Believe: at n = 7 class 2, 8 of 9 E-151 failing steps have parent ~ child (Hom-verified F/R chains to seed LNAs; the two seeds are equal, isomorphic, or joined by <= 5 edges). Record 6 has no child path. Class 1 and the other 28 class-2 E-204 edges not done. `073/toolsmith_certificate.py` (modes chain, corrupt, edges, seedlinkiso, f010), needs `049/toolsmith_collect.py 7 2 20000 pkl 100` (311 s).
- Believe: `canonicalKey` is a LABELLED key: isomorphic seeds (04000 and its relation dual, 4 more pairs in class 2) have different keys, so labelled meets can miss. Parent and child of a step may reach different seeds; E-157's "joined to an LNA" did not say which.
- Believe: no NO-test has demonstrated power at n = 7: key-equal inequivalent pairs do not exist there; cross-key rejection is forced by the Coxeter polynomial. Hom test: 0 disagreements over ~150k verdicts, key-moving 52/52 rejected (shards overlap). Cartan-level invariants (SNF of C, C+C^T, C-C^T, signature) equal for the F-010 pair, so no power there either.
- Round 073 review: NO-power wording narrowed (1 Hom rejection in 195 corruptions; rest legality/iso). `parallel-at-meet` False at all 16 meets. Overnight script `073/toolsmith_f010_overnight.py --budget-hours`.
- Live NO control: n = 9 F-010 with iso meets: 4 + 4 no meet (116 s), 5 + 5 YES control 507 s, NO control exceeds 10 min (OVERNIGHT proposal).
- Tried and dropped: n = 7 key-moving steps from tilting-closed neighbourhoods (0 of ~800); they only occur deep in the ball.
- Next: collect the remaining 28 class-2 and 34 class-1 edges (collector caps at 20 000 expansions; a depth-8 complete collector is 8x slower); make iso-meet the default in `meetingPoints`; F-010 5 + 5 overnight.
- Traps: shards of `keymoving` repeat levels 1-7, so per-shard tallies must not be summed; shell loops over many runs hit the 120 s foreground limit and go to background (use files).
