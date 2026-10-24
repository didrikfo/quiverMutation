# Toolsmith notebook

## Round 003-065 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED; meet up to sign-aware iso. End(T) iso at every step: only Hom(T,T[-1]) tests the premise.

## Round 067-081 (T10/T5)
- E-168/E-198: 70 key-kept failing edges; 073-078: n = 7 failing steps joined by Hom-verified chains (`049/toolsmith_collect.py`, `078/toolsmith_chains.py`); `CyclicDescendantsError`/`allowCyclic` (E-211).
- 081: class 2 seeds 1 and 8 joined (`081/toolsmith_seedlink.py`); use a ball per seed. Hom+Cartan has no NO power (E-228); the iso meet is the check. Trap: exec'd preambles use global `P`.

## Round 083
- 3060000 ball depth 5: 10 key-kept perI J != 0 steps (10/4770), all Hom-real, all children join the root at 7-8 (E-234).

## Round 085
- Believe: 3304000 ball depth 5 (4 046 nodes, 2 958 parents): 1 key-kept perI step (index 2026, F2 R8 R9 R7 R7, R7; Hom(T,T[-1]) = 2, fails tiltingPlus), child joins 3304000 at total 6. Rate 0.034% vs 0.21% for 3060000 (one event, no test). Both sides: key-kept J != 0 steps exist at depth 5, none leaves its own J = 0 class.
- The 10 E-234 children are 8 iso classes ({1915,4722,E-233 child}, {1970,2003}); 3304000's child iso to none of the 11 (consistent with F-010 inequivalence, a free sanity check).
- Tools: `085/toolsmith_d5scan.py` now takes env KEYROOT, DEPTH (key base was hard-coded); keys of 3060000 and 3304000 asserted equal. Builds: 3304000 depth 5 116 s, 3060000 131-170 s; scan ~35 s per 1500 parents. `toolsmith_distinct.py` = pairwise iso of children from `d5pkl` pickles.
- Traps: `pgrep -f` in a monitor matches its own shell, loops never end; sandbox wall-clock and my sense of time disagree, poll with a bash `until` under 115 s.
- Next: depth-6 scan (3304000 build ~8 min, borderline; 3060000 needs checkpointing: propose for OVERNIGHT); skeptic recheck of index 2026; is a J != 0 step ever followed by a second one (J != 0 child as parent)? a NO-capable certificate remains the real gap, these steps never test the premise since their joins are J = 0.
