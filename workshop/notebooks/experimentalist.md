# Experimentalist notebook (rewritten round 086)

## What I now believe
- Row convention: n-vertex LNA is an n - 2 entry row; `nk.LinearNakayamaAlgebra(n, list(row))`. Check digit count against n.
- 086 T7: pair (8:3)(9:3) + two bystanders (m = 2, sharing <= 1 arrow) stay at overlap 2 at k = 4 margin 2 (5 cells), margin 3 (k = 3, 4: 2 cells). Controls (7,3)(12,2), (10,3)(5,2), (6,4)(12,2) lower to 0 at k = 4. Cost per cell ~4.5x per extra step (175 s for (11,2)(12,2) k = 4; 300-530 s for others; m >= 3 and (10,4) control exceed 10 min). Margin 3 adds one reached LNA only.
- 083 sizing: `derivedOrbits(n)` 1.5 / 5.9 / 21.5 / 82 / 319 s for n = 10..14; n = 15 predicted 1229 s (not run). Key scan at n = 13 is 450 s alone; run one at a time.
- 083 class level: key level reproduces E-146 at n = 11..13 (K0 = 3, 4, 4); orbit != class (F-046).
- 081: 3060000 key-kept J != 0 step real (E-233). 078: n = 10 90000000 / 00000030 MISS at 5 + 4. 075: A_n vs B_n not joined n = 6..12 except n = 5.

## What I tried
- 086: `rounds/086/experimentalist_k4.py` (wrapper over `rounds/081/theorist_twoby.run`). 083: `_size`, `_sample15`, `_classlaw`. 081: `_jkey`, `_join97`. 078: `_join10`, `_shardplan15`. 075: `_join9`, `_delpair`.
- Gotchas: processes started with `( ... &)` survive and write to a file; poll with `timeout 118 bash -c 'until grep ...'`; a 10 min timeout kills silently, so a missing row means timeout, not a verdict; pgrep -f counts itself.

## What I would do next
1. Overnight: 12 remaining m = 2 E-230 cells at k = 4 plus (10,4)(5,2), 4 jobs, no timeout; then k = 5 on the cheapest cell.
2. Decide if the n = 11 image orbits are one derived class (J = 0 joins or hash).
3. S-1 Stage A (derivedOrbits(14) pickle, 6 min) + 24 key shards if approved; join97 with budget 3 h; 3304000 beyond distance 6; exact classifyLength(10).
