# Experimentalist notebook (rewritten round 091)

## What I now believe
- Row convention: n-vertex LNA is an n - 2 entry row; `nk.LinearNakayamaAlgebra(n, list(row))`. Check digit count against n.
- 091 (T10 calibration): 2400230 ball (E-212 moves) sizes 1, 9, 40, 130, 369, 1001, 2691, 7261, 20210 (depth 0..8; 3304000 is ~4x larger at depth 5). 1 key-kept perI step in 20210 parents (depth 8, child joins its own ball at 9). (Flags only from depth 6: 3, 47, 348; key(2400230)=key(3060000)=key(3304000)=base.) No bridge to 3304000 among explored radii: C(4) vs R(5) none, A(6) vs R(5) none. Recorded J != 0 children lie in their roots' J = 0 components, so adding them cannot bridge; only 8427 and unscanned exits could.
- 089 (S-1): |orbit(B_n)| = n 2^(n-5) for n = 6..16; n = 17 predicted 69632. Job 25 s via `rounds/089/experimentalist_orb13.py`.
- 086 T7: (8:3)(9:3) + 2 bystanders stay at overlap 2 at k = 4; cost ~4.5x per step (175-530 s per cell).
- 083: `derivedOrbits(n)` 1.5 / 5.9 / 21.5 / 82 / 319 s n = 10..14; key scan n = 13 450 s.
- 081: 3060000 J != 0 step real (E-233). 078: n = 10 90000000 / 00000030 MISS at 5 + 4.

## What I tried
- 091: `jscan` (scan whole ball for kept perI), joins via 082/skeptic_join.py side/meet. 089: `orb13`. 086: `k4`. 083: `_size`, `_classlaw`. 081: `_jkey`, `_join97`. 078: `_join10`. 075: `_join9`.
- Gotchas: no pickles of earlier balls are in the repo (they live in /tmp), so every ball is rebuilt; `build` with DEPTH env and `--budget-hours 0.14` exits 0 when budget spent (read the log, rerun); piping through tail hides progress, use a file; the auto-mode classifier sometimes fails transiently, retry once; a 10 min timeout kills silently, so a missing row is not a verdict; pgrep -f counts itself.

## What I would do next
1. Overnight: 2400230 ball depth 9-10 + jscan; 8427 join; meet 2400230(8) vs 3304000(6).
2. Other relation placements of B_n (orbit counts are cheap).
3. S-1 Stage A (derivedOrbits(14) pickle, 6 min) + 24 key shards if approved; 12 m = 2 E-230 cells at k = 4.
