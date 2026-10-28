# Experimentalist notebook (rewritten round 094)

## What I now believe
- Row convention: n-vertex LNA is an n - 2 entry row; `nk.LinearNakayamaAlgebra(n, list(row))`. Check digit count against n.
- 094 (T7): k = 5, margin 1, pair (8:3)(9:3): 5 of 5 m = 2 cells stay at overlap 2, controls (7,3)(12,2), (10,3)(5,2) drop to 0. Margin 2 at k = 5 does not fit 10 min (est 790-2400 s per cell, x4.5 per step). Reached counts at k = 5 margin 1: 31..57.
- 091 (T10): 2400230 ball (E-212 moves) sizes 1, 9, 40, 130, 369, 1001, 2691, 7261, 20210 (depth 0..8). No bridge to 3304000 among explored radii; key(2400230)=key(3060000)=key(3304000)=base.
- 089 (S-1): |orbit(B_n)| = n 2^(n-5) for n = 6..16; n = 17 predicted 69632 (25 s via `rounds/089/experimentalist_orb13.py`).
- 086 T7: k = 4 margin 2 cells 175-530 s each. 083: `derivedOrbits(n)` 1.5/5.9/21.5/82/319 s n = 10..14. 078: n = 10 90000000 / 00000030 MISS at 5 + 4.

## What I tried
- 094: `k5` (plan wrapper), `k5par` (first-move split of localMutationSequences; no real speed-up, 207 s vs 175 s at k = 4, because the shared `seen` table is lost). Margin 1 was the way to fit k = 5.
- Earlier: `jscan`, `orb13`, `k4`, `_size`, `_classlaw`, `_jkey`, `_join97`, `_join10`, `_join9`.
- Gotchas: no pickles of earlier balls in the repo (they live in /tmp); piping through tail hides progress, use a file; a 10 min timeout kills silently (exit 124), so a missing row is not a verdict; a foreground command over 600 s is moved to background by the harness, the timeout still kills it; pgrep -f counts itself; two cells at 2 jobs each is the efficient packing on 4 cores.

## What I would do next
1. Overnight: T7 k = 5 margin 2 (8 cells serial, ~1 h on 4 cores); 2400230 ball depth 9-10 + jscan; 8427 join.
2. Run margin 1 at k = 4 to see whether margin 1 loses anything vs E-241.
3. Other pair shapes at k = 3 (E-152 sampled three); S-1 Stage A if approved.
