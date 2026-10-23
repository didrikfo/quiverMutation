# Experimentalist notebook (rewritten round 083)

## What I now believe
- Row convention: n-vertex LNA is an n - 2 entry row; `nk.LinearNakayamaAlgebra(n, list(row))`. Check digit count against n.
- 083 S-1 sizing: `derivedOrbits(n)` (default free=True) measured 1.5 / 5.9 / 21.5 / 82 / 319 s for n = 10..14; rows 4862 / 16796 / 58786 / 208012 / 742900; orbits 746 / 3712 / 14274 / 52257 / 188093; n = 14 peak RSS 467 MB. Sampling the loop body on 3000 rows predicts each within 6%; n = 15 = 1229 s (not run).
- 083 class-level law: key level reproduces E-146 at n = 11, 12, 13 (K0 = 3, 4, 4). Orbit labels of images split (2-3 orbits at the K0 level, 1 key), ext orbits (edges+doubles) 2 at n = 11 (2,5) and n = 13 (4,5)/(2,7). Undecided whether those are one class. Orbit level cannot confirm or refute the law.
- Key scan at n = 13 is 450 s alone (2.2 ms/row); n = 13 `classlaw --ext` takes 527 s of the 10 min. Do not run more than one at a time.
- 081: 3060000 had one key-kept J != 0 step (E-233 settled it: real, child joined at 7). 078: n = 10 90000000 / 00000030 MISS at 5 + 4; n = 15 key plan 24 shards x 3.8 min after an n = 14 table.
- 075: A_n vs B_n family n = 6..12 not joined by orbit + mirror except n = 5.

## What I tried
- 083: `rounds/083/experimentalist_size.py`, `_sample15.py`, `_classlaw.py`. Gotchas: Bash blocks `sleep N; cat`; wall clock runs slower than Monitor timers suggest, use `timeout 115 bash -c 'until grep ...'`; ToolSearch select Monitor before use; processes started with `( ... &)` survive.
- 081: `_jkey.py`, `_join97.py`. 078: `_join10.py`, `_shardplan15.py`. 075: `_join9.py`, `_delpair.py`.
- Lessons: grep EXPERIMENTS.md first; mirror needs len(row)+2; forward orbit absent != not equivalent; tee to a file.

## What I would do next
1. Decide if the n = 11 image orbits (3 at (3,4) K >= 3 -> n = 10; 2 at (2,5)) are one derived class (cheapest case), by J = 0 joins or a hash.
2. If approved: S-1 Stage A (derivedOrbits(14) pickle, 6 min) + 24 key shards, then Stage C.
3. Overnight: join97 with --budget-hours 3; 3304000 beyond distance 6; n = 10 A-ball pickle.
4. Exact classifyLength(10) overnight.
