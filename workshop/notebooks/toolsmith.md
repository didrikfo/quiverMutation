# Toolsmith notebook

## Rounds 003-102 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s (use `timeout 115 sh -c 'while pgrep ...'` loops), `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); Hom+Cartan has no NO power for iso, the iso meet is the check.
- T10 n = 9 (E-234..E-268): key-kept J != 0 steps exist; join searches cannot test the J = 0 premise. Tools: `087/toolsmith_d6scan.py`, `098/toolsmith_j98.py`, `102/toolsmith_litex.py`.

## Rounds 105-117 (n = 7 class 1; n = 8, 9 q_S)
- "A_n" = Nakayama rad^(n-1)=0, "B_n" = lone 3-arrow relation; both ~ kD_n. All 13 key-kept steps of E-201 Hom-verified (E-297, E-300); meet iso does the rejecting. Tool `114/toolsmith_s8par.py`.
- R117: n = 8 no key-equal inequivalent LNA pair; n = 9 one split group (F-010); q_S is NO-capable but cannot fire on J = 0 steps (derived equivalences keep q_S), T10 premise untested. `117/toolsmith_qsscan.py`, `toolsmith_qsball.py` (needs `DEPTH`; sets `base` itself).
- Traps: `moves`/`jcands` use global `base`; /tmp checkpoints vanish between sessions.

## Round 118 (backlog 11, 12: both settled)
- Docstrings edited (search guard paragraph, canonicalKey + fingerprint module, mutateAtVertex allowCyclic + Oppermann note); 60 tests pass. F-002 rewording proposed to the chair, not edited.
- Fosse Cor. 4.6 misprint re-derived from LaTeX (only v1 exists): rows 3, 4 should end n_0 + 2. = E-285 rediscovery.
- Five direct-caller scripts: three exact, 034 and 027 are wall-clock capped; under load they cover less (034: 7 435-7 805 vs 8 456 expansions; 027: 4 445-6 333 vs 9 197) but agree on the covered part. Trap: do not compare capped-walk counts across runs without checking load; run alone, and never four jobs at once.

## Would do next
- Run qsball at depth 6 for 3060000 (checkpointed, overnight): does a key-kept step cross to the 3304000 q_S half? (T10 premise.)
- n = 10 key groups with q_S splits (E-292 has 2): sizing.
- If asked: uncapped rerun of 034 tail (3 ROW lines, ~12 min alone) with a larger budget in two runs.
