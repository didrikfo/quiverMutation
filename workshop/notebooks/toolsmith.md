# Toolsmith notebook

## Round 003-065 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED; meet up to sign-aware iso. End(T) iso at every step: only Hom(T,T[-1]) tests the premise.

## Round 067-081 (T10/T5)
- E-168/E-198: 70 key-kept failing edges; 073-078: n = 7 failing steps joined by Hom-verified chains (`049/toolsmith_collect.py`, `078/toolsmith_chains.py`); `CyclicDescendantsError`/`allowCyclic` (E-211).
- 081: class 2 seeds 1 and 8 joined (`081/toolsmith_seedlink.py`); use a ball per seed. Hom+Cartan has no NO power (E-228); the iso meet is the check. Trap: exec'd preambles use global `P`.

## Round 083
- Believe: in the 3060000 ball (n = 9), depth-5 parents hold 10 key-kept perI J != 0 steps (10/4770 = 0.21%, 10 of 640 perI candidates), all fail tiltingPlus with Hom(T,T[-1]) = 2-4, all children join the root at 7-8 by J = 0 moves with iso meet. Plus E-233's depth-4 one = 11. 3304000 not scanned at depth 5.
- Tools: `083/toolsmith_d5scan.py build/scan/ctl` (ball 167 s, scan of 4770 parents ~110 s total: the scan is cheap, the ball build is the cost), `toolsmith_d5pkl.py` makes skeptic_kept-format pickles so `082/skeptic_kept.py test` and `082/skeptic_join.py` run unchanged.
- Not done: 3304000 depth-5 scan; depth-6 parents (~19k nodes, ball build > 10 min, needs checkpointing); distinctness of the 10 children; meet totals were 8 for six (child depth 3 only; deeper child ball may give 7); first five join outputs only in the terminal.
- Next: 3304000 scan (same scripts, root arg and key base need parameterising: `base` is hard-coded to 3060000's key, same for both); checkpointed build for depth 6; a negative control for the verifier on these joins.
