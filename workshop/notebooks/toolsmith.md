# Toolsmith notebook

## Rounds 003-098 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); meet up to sign-aware iso. Hom+Cartan has no NO power, the iso meet is the check.
- T10 n = 9 (E-234..E-255, round 098): key-kept J != 0 steps exist; 2400230 vs 3304000 calibration exhausted at feasible radii (balls share no WL hash), a move-set statement, nothing on the J = 0 premise. Not an OVERNIGHT case.
- Tools: `087/toolsmith_d6scan.py`, `090/toolsmith_nonekey.py`, `093/toolsmith_unionjoin.py`, `098/toolsmith_j98.py`. Traps: scripts exec the 087 prelude; `background &` makes task notice fire at once, poll logs.

## Round 102 (literature, T5/T10)
- Ran the worked examples of 2112.08129 (Ex. 3.1, 3.2, Sect. 6) and 1504.02617 (Ex. 4.1, 4.2) through gate / J (`perI`, copied from rounds/033) / key guard: `rounds/102/toolsmith_litex.py` (32 s), output `toolsmith_litex.out`, submission `rounds/102/toolsmith.md`.
- Agree: Fosse Ex. 3.1, Opp Ex. 4.1, Opp Ex. 4.2 (repo result = degree-0 part of the printed final, iso to A^op); Sect. 6 n = 3, 4 reach the 2 x n rectangle in 5 / 12 gate-admitted, J = 0, key-kept steps (BFS, not the paper's sequence). Disagree: Fosse Ex. 3.2 (cyclic): refused (E-211), and `allowCyclic=True` is a `KeyError` (procedure.py l. 221), not a wrong quiver.
- Oppermann's LEFT mutation = repo RIGHT mutation in these labels; repo `left` is refused by gate and key guard (a free power control).
- Oppermann v1 prints no loop-of-degree->=1 example (only Ginzburg loops, theorem level); library cannot represent degrees.
- Mine, unverified: summary's "Fosse step 3 cyclic <-> Opp step 4" is really Opp Step 2 with phi the degree-1 loop (+ Step 4 adds alpha omega alpha^-1).
- Traps: `perI` on a cyclic quiver returns garbage (`{2: -1}`); `bfs` over F moves with key guard reaches n = 4 rectangle at depth 9 in 32 s.

## Next
- Ex. 3.2 by the dg rule by hand (theorist); `allowCyclic=True` clean error (small, with test) only if asked.
- 3304000 / 2400230 deeper scans only widen a MISS; do them only if asked.
