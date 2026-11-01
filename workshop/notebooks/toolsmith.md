# Toolsmith notebook

## Rounds 003-102 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); meet up to sign-aware iso. Hom+Cartan has no NO power for iso, the iso meet is the check.
- T10 n = 9 (E-234..E-268): key-kept J != 0 steps exist; join searches cannot test the J = 0 premise.
- Round 102: worked examples of 2112.08129 / 1504.02617 through gate/J/key guard (`rounds/102/toolsmith_litex.py`); Fosse Ex. 3.2 cyclic refused, `allowCyclic=True` is a KeyError. Oppermann LEFT = repo RIGHT.
- Tools: `087/toolsmith_d6scan.py`, `098/toolsmith_j98.py`, `102/toolsmith_litex.py`. Scripts exec preludes; `background &` fires task notice at once.

## Round 105 (H-027 / T10, independent tilting test)
- The independent test already exists: `rounds/050/skeptic_tilt.py` (`Alg`, `stepTest(alg,k,S)`, S=-1), E-161. Reuse it (exec it); do not rewrite.
- `rounds/105/toolsmith_tiltA6.py`: all 10 one-step complexes of B_6 (right, left = opposite algebra): 8 tilting, 2 not, matching the gate exactly; no H equals Cartan(A_6) up to perm/transpose. So A_6 is not a one-step mutation of B_6 (necessary-condition argument). Gate agreement 281/281 on all LNA n=4..6; planted wrong shift and repeated summand rejected (det H = 0).
- `toolsmith_tilt8.py`: 7 of 8 E-268 J != 0 edges (n=9, replayed from printed paths, no ball needed) all have Hom(T,T[-1]) != 0, H = Cartan(child)^T. Known pattern (E-161/E-184); the 8th (E-239, depth-5 B scan) not rebuilt.
- Believe: J != 0 steps are silting-not-tilting; Cartan(End T) = child's Cartan holds anyway, so Cartan equality is no evidence of equivalence.
- Trap: left step = right step on `a.opposite()` (labels kept, vertex 1 has no out-arrows there); the first run printed the gate for the wrong side (fixed). `repoChild` must use `isMutable`, not `mutateAtVertex` (does not gate).

## Next
- If A_6 from B_6 is wanted: Cartan-first search for a tilting complex (summand = 2-term complexes) with H ~ Cartan(A_6), at n = 6; Fosse Cor. 4.6 for the explicit T.
- 8th J != 0 edge (E-239) needs the 3304000 depth-5 scan (~2 min after a ~11 min ball build); only if asked.
- 3304000 / 2400230 deeper scans only widen a MISS.
