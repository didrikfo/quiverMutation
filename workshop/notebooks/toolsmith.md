# Toolsmith notebook

## Rounds 003-102 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell (happened again r107), > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); meet up to sign-aware iso. Hom+Cartan has no NO power for iso, the iso meet is the check.
- T10 n = 9 (E-234..E-268): key-kept J != 0 steps exist; join searches cannot test the J = 0 premise.
- Round 102: worked examples through gate/J/key guard (`rounds/102/toolsmith_litex.py`); Fosse Ex. 3.2 cyclic refused. Oppermann LEFT = repo RIGHT.
- Tools: `087/toolsmith_d6scan.py`, `098/toolsmith_j98.py`, `102/toolsmith_litex.py`. Scripts exec preludes; `background &` fires task notice at once.

## Round 105 (E-282)
- `rounds/050/skeptic_tilt.py` (`Alg`, `homdim`, `stepTest`) is the independent chain-map Hom test; `105/toolsmith_tiltA6.py` = one-step B_6 test (8 tilting, 2 silting-only, no A_6 Cartan). J != 0 n=9 edges are silting, not tilting.
- Trap: left step = right step on `a.opposite()`; `repoChild` must use `isMutable`.

## Round 107 (H-027 multi-step)
- "A_n" = Nakayama rad^(n-1)=0, "B_n" = lone 3-arrow relation; both derived equivalent to kD_n (Euler form has n(n-1) roots: 20, 30), NOT type A.
- All indecomposables of K^b(proj B_n) up to shift = n(n-1); brute-force thin complexes miss 7 / 12; closing under cones of chain maps (`107/toolsmith_cones.py`, minimise by Gaussian elimination) finds all in seconds. Exceptional test (End=k, no self-ext) = indecomposable.
- `107/toolsmith_tiltsearch.py`: Hom table h[r][s][j] once, then tilting-exchange BFS is pure table lookup (n=6 depth 5: 1745 nodes, 15 s). Result: A_n reached at depth 3 (n=5) / 4 (n=6), End T = A_n verified by composites of arrow maps (homotopy-class test with controls). One-step: 8 vertex-local (= E-282) + 2 non-local at n=6.
- Bug caught: a coupled |shift| filter in compat() gave 6 complements for an almost complete T (should be 2); the count of complements per summand is a good sanity check. Believe: generation = rigid + n summands + det +-1 (Bruning, UNVERIFIED).
- Belief: E-280's "A_n not in orbit(B_n)" is a statement about gate-limited LNA moves, not about tilting reachability.

## Next
- n = 7, 8, 9 same machinery (42..72 objects); depth-to-first-hit vs n; shifts window check.
- Describe End of each intermediate T on the depth-4 path (LNA? which step does the gate refuse?) with a theorist; compare with Fosse's explicit T.
- Verify exchange = Aihara-Iyama mutation (silting order) before it enters a finding.
