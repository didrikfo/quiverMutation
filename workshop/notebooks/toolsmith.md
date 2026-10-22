# Toolsmith notebook

## Round 003-065 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s (`timeout N bash -c 'until ...'`), `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (None for bundles > 720 relabelings); meet up to sign-aware iso. End(T) iso at every step: only Hom(T,T[-1]) tests the premise.

## Round 067-078 (T10/T5)
- E-168/E-198: 70 key-kept failing edges. 073: n = 7 class 2, 8 of 9 E-151 failing steps joined by Hom-verified chains (`049/toolsmith_collect.py`, 311-405 s). 075: iso-meet joins of the 4 unjoined E-213 classes. 078: all 8 chains recorded and Hom-verified (`078/toolsmith_chains.py R 4 5`); `CyclicDescendantsError`/`allowCyclic` in `procedure.mutateAtVertex` (E-211).

## Round 081
- Believe: seeds 1 and 8 of class 2 (n = 7) are joined by F6 | R5 R5 R6 (4 edges, Hom-verified, iso meet); closes E-222 rep 2 and E-212 record 6 (child -> seed 1, parent -> seed 8, same pair). Pooled target ball hides seed-to-seed links; use a ball per seed (`081/toolsmith_seedlink.py 2 3 3`, seconds).
- `081/toolsmith_rec6.py PKL 6 4 6` = chains script on a pickled record (~100 s after the 405 s collect). Trap: exec'd preambles use global `P` as a prime; never name algebras `P`/`C` there.
- Not done: composite chain replayed end to end; class 1 records; negative control for the verifier (still thin, E-212); all 14 seeds pairwise.
- Next: pairwise seed-link matrix for classes 1 and 2 (cheap); F-010 5 + 5 overnight proposal; iso-meet default in `meetingPoints`; DEFAULT_CAP rewording.
