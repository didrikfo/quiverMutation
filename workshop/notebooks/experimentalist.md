# Experimentalist notebook (rewritten round 089)

## What I now believe
- Row convention: n-vertex LNA is an n - 2 entry row; `nk.LinearNakayamaAlgebra(n, list(row))`. Check digit count against n.
- 089 (S-1): |orbit(B_n)| = n 2^(n-5) holds n = 6..16 (3328 at 13, 7168 at 14, 15360 at 15, 32768 at 16 in 34 s), closed, mirror-closed, A_n not in it. Whole job 25 s via `rounds/089/experimentalist_orb13.py`. Fit still unproved; n = 17 predicted 69632.
- 086 T7: pair (8:3)(9:3) + two bystanders (m = 2) stay at overlap 2 at k = 4 margin 2; controls lower to 0. Cost ~4.5x per extra step (175-530 s per cell; m >= 3 and (10,4) exceed 10 min).
- 083 sizing: `derivedOrbits(n)` 1.5 / 5.9 / 21.5 / 82 / 319 s for n = 10..14; n = 15 predicted 1229 s (not run). Key scan at n = 13 is 450 s alone.
- 083 class level: key level reproduces E-146 at n = 11..13; orbit != class (F-046).
- 081: 3060000 key-kept J != 0 step real (E-233). 078: n = 10 90000000 / 00000030 MISS at 5 + 4. 075: A_n vs B_n not joined n = 6..12 (now ..16 by orbit).

## What I tried
- 089: `orb13`. 086: `k4`. 083: `_size`, `_sample15`, `_classlaw`. 081: `_jkey`, `_join97`. 078: `_join10`, `_shardplan15`. 075: `_join9`, `_delpair`.
- Gotchas: `( ... &)` processes survive and write to a file; poll with `timeout 118 bash -c 'until grep ...'`; a 10 min timeout kills silently, so a missing row is a timeout, not a verdict; pgrep -f counts itself. orbitReport is fast (seconds at n = 15); orbit counts are cheap, derived-class work is the expensive part.

## What I would do next
1. Overnight: 12 remaining m = 2 E-230 cells at k = 4 plus (10,4)(5,2), 4 jobs, no timeout; then k = 5 on the cheapest cell.
2. Other relation placements of B_n (3-arrow relation not at the end): orbit counts are cheap, test the formula's dependence.
3. Decide if the n = 11 image orbits are one derived class; S-1 Stage A (derivedOrbits(14) pickle, 6 min) + 24 key shards if approved.
