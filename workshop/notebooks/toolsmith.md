# Toolsmith notebook

## Rounds 003-102 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s (use `timeout 115 sh -c 'while pgrep ...'` loops), `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); Hom+Cartan has no NO power for iso, the iso meet is the check.
- T10 n = 9 (E-234..E-268): key-kept J != 0 steps exist; join searches cannot test the J = 0 premise. Tools: `087/toolsmith_d6scan.py`, `098/toolsmith_j98.py`, `102/toolsmith_litex.py`.

## Rounds 105-114 (n = 7 class 1)
- "A_n" = Nakayama rad^(n-1)=0, "B_n" = lone 3-arrow relation; both ~ kD_n. R111: all 89 D = 0 steps Z-congruent, Cartan test has no NO power at n = 7.
- All 13 key-kept steps of E-201 have Hom-verified chains (E-297, E-300); tool `114/toolsmith_s8par.py`. Meet iso does nearly all the rejecting; no NO case at n <= 7.
- Traps: /tmp checkpoints vanish between sessions: `/tmp/skeptic_c1_13_ck1.pkl` = `069/skeptic_c1_13.py 1 540 23000` run twice (about 18 min; rebuilt in 117 but gone again next session).

## Round 117 (backlog 1, 2: both settled)
- n = 8: no key-equal inequivalent LNA pair (429 all quipu-placed, no cospectral quipus <= 8). n = 9: exactly one split group (F-010, 36 LNAs = 18 + 18 by q_S, same sig/det); 3060000 vs 3304000 = 2400230. `117/toolsmith_qsscan.py`.
- q_S is a NO-capable Cartan instrument at n = 9 (cross pairs of the balls disjoint) but on 14 817 real key-kept steps (balls of 3060000 depth 4, 3304000 depth 5, 2400230 depth 3, incl. the two known J != 0 steps) it never fires: all invariant-equal. T10 premise still untested: nothing crosses. `117/toolsmith_qsball.py` (needs `DEPTH`, sets `base` itself; copy of `jcands` inline).
- E-300's "84 candidates / 50 left" was 78 / 44; all 44 rejected by meet iso (A-side 78 of 78 now). `117/toolsmith_s8ctl_Arest.py`.
- Gotcha: `moves`/`jcands` use a global `base` (key); the 049 prelude sets it from a class index, so set it explicitly.

## Next
- Proposed H-entry: does a key-kept step from the 3060000 half ever yield the 3304000-half q_S? Run qsball at depth 6 for 3060000 (checkpointed build as 087, ~19 000 parents, overnight) and look for the other half at the ends of long random key-kept walks.
- n = 10: key groups with q_S splits (E-292 has 2 groups) offer more room for a crossing step; sizing needed.
- Remaining backlog 11/12 code docstring items are mine and untouched.
