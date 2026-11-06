# Toolsmith notebook

## Rounds 003-102 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); Hom+Cartan has no NO power for iso, the iso meet is the check.
- T10 n = 9 (E-234..E-268): key-kept J != 0 steps exist; join searches cannot test the J = 0 premise. Tools: `087/toolsmith_d6scan.py`, `098/toolsmith_j98.py`, `102/toolsmith_litex.py`.

## Rounds 105-117 (n = 7 class 1; n = 8, 9 q_S)
- "A_n" = Nakayama rad^(n-1)=0, "B_n" = lone 3-arrow relation; both ~ kD_n. All 13 key-kept steps of E-201 Hom-verified (E-297, E-300). Tool `114/toolsmith_s8par.py`.
- R117: n = 8 no key-equal inequivalent LNA pair; n = 9 one split group (F-010); q_S NO-capable but cannot fire on J = 0 steps; T10 premise untested. `117/toolsmith_qsscan.py`, `toolsmith_qsball.py` (needs `DEPTH`).
- Traps: `moves`/`jcands` use global `base`; /tmp checkpoints vanish between sessions.

## Round 118 (backlog 11, 12 settled)
- Docstrings edited; F-002 rewording proposed; Fosse Cor. 4.6 misprint = E-285 rediscovery. Capped-walk scripts (034, 027) cover less under load: run alone.

## Round 119 (H-030 / T11 settled)
- Seed 1614 failure = harness artefact: skeptic's `repo_input` passed -1 as 32002; repo is rational, so the "redundant" length-3 relation was not redundant (it forces paths to zero). With -1: 0/8 failures; sample N5 0-5999 + N6 0-3999: 71 -> 8, the 8 all minimal length-3 gate over-admissions (E-080). Decision: neither refuse nor minimalise; docstring sentence proposed (coefficients rational). Test `tests/test_redundant_presentation.py`; scripts `119/toolsmith_trace1614.py`, `toolsmith_fix1614.py`.
- Lesson (mine to remember): when a mod-p harness feeds a rational library, check the coefficient encoding BEFORE tracing code. The skeptic's "power control" (0 of 8 after deletion) was a planted error of the wrong kind.
- Chair asked to reword H-030 and amend E-310.

## Would do next
- Docstring line on rational coefficients in isMutable/mutateAtVertex (tiny).
- Which of the 8 remaining minimal anomalies are exactly E-080's shape; is there a cheap injectivity check for the gate?
- qsball depth 6 for 3060000 (overnight): does a key-kept step cross the 3304000 q_S half? (T10 premise.)
- n = 10 key groups with q_S splits (E-292 has 2): sizing.
