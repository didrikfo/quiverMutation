# Toolsmith notebook

## Rounds 003-102 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); Hom+Cartan has no NO power for iso, the iso meet is the check.

## Rounds 105-119 (short)
- n = 7..9 key-kept steps Hom-verified (E-297, E-300); q_S instrument (`117/toolsmith_qsscan.py`, `qsball.py`). R119: seed 1614 = harness coefficient encoding (-1 as 32002), E-312. Lesson: check encodings before tracing code.

## Round 122 (T10): FIRST CROSSING FOUND (pending review)
- Tool `122/toolsmith_jxsearch.py`: `--plan N` lists INV-split key groups (INV = SNF(C+C^T), signature, q_S hist), cached in `toolsmith_jxgroups.json` (n = 9: 1, n = 10: 3, n = 11: 6 groups; more than E-292 because it includes singular S and signature splits). `run N` does a guarded F/R BFS that expands every key-kept child, labels J by `perI`, and union-finds J = 0 components (up to iso via 075's WL hash + explicit iso, ISO=1). For each key-kept J != 0 step it records INV(parent) vs INV(child), joined/OPEN, and an EXT-radius J = 0 extension. CTRL=planted grafts a crossing edge and an in-component edge.
- Hit: n = 11 group 5 (72 = 13 + 59, SNF tails 2,4 / 2,2), root 003033030. P = F3 R8 F3 (J = 0); the F3 step at P is J != 0, key-kept, gate-admitted, and SNF goes 2,4 -> 2,2 (det C = 1, acyclic). The child is one J = 0 step (F6) from LNA 030330030 (59-class), so a 5-step guarded F/R walk joins two derived-inequivalent LNAs: H-015 false. Checked with brute-force Cartan and determinantal divisors in `toolsmith_jxhit.py`. The library's F-only walk (`toolsmith_jxlib.py`) does not reach it at depth <= 6.
- Null elsewhere: 19 other roots (n = 9 depth 5, n = 10, 11 depth 4) have 0 key-kept J != 0 steps, except E-233's pair at 3060000, which joins up to iso at ext 2. Labelled keys made it look OPEN: always use ISO.
- No flint/pari in .venv: use sympy `smith_normal_form` (fast enough at n <= 11). The n = 11 group table costs about 5 min: use the cache.
- Belief now: the key guard is necessary only; the J != 0 refusal (E-151, 8-11% cost) should become a library guard. Merges resting on J != 0 steps need re-audit.

## Would do next
- After review: library opt-in `refuseNonTilting=True` in `mutationSearchDepthFirst` (J check = perI = tiltingPlus), with a test built on P -> X from 003033030 (gate admits, key kept, refused under the option). Time it.
- Audit: which recorded merges (E-032, F-entries at n >= 9) used a J != 0 step? Scan replayable paths with the perI flag.
- n = 11 group 5: all 72 roots at depth 4; F-only depth 7-8 from 003033030 (does the library walk itself cross?).
- n = 12 sizing (group table > 5 min; propose overnight). Docstring line on rational coefficients (still pending). E-080 shape of the 8 anomalies (T11).
