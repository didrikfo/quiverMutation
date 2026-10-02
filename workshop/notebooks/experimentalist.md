# Experimentalist notebook (rewritten each round)

## What I now believe (after round 025)
- 025 (capped 500 s walks, `rounds/023/scholar_longsquare.py`, 4 parallel runs, files `rounds/025/experimentalist_n*_c*.txt`): J != 0 with out-degree 2 and no long square: n = 8 c0 48 (10 478 algebras), n = 8 c1 15 (16 976), n = 8 c3 0 (18 269, only 6 J != 0 rows), n = 9 c0 0 (8 231 algebras, 85 J != 0 rows all out-1 + long square). Counts are cap lower bounds; zeros are weak.
- New oddity: n = 8 c1 has 2 J != 0 rows with out-degree 1 and NO long square (not inspected; maybe presentation artefact).
- J != 0 <=> tiltingPlus False held on every row of all runs. Only the long-square shape breaks (n = 8).
- T5 (022): n = 5..7 long square at 100 % of rejects, 0 of 479 761 tilting steps. Holds at n <= 7 only.
- T5 (019): dim ker 0 on tilting steps, >= 1 on non-tilting. From 017: walk tables robust to the reduceAgainstPivots fix; `MAXEXP=N` caps deterministically. From 014: n = 17 orbits closed.

## What I tried
- 025: four parallel `--budget-sec 500` runs; parallel runs cover fewer algebras than solo (referee 42 in 200 s solo vs 48 here in 500 s parallel). Prefer solo or deterministic caps.
- 022: shape tests on every step; wall-clock caps make counts differ between runs: quote rates.

## What I would do next
1. Print one example per tally key and classify the n = 8 c1 out-1 no-longsq rejects and the out-2 rejects (D vs G).
2. Solo longer runs: n = 9 c0, c1 with 3000 s overnight; n = 8 c2; deterministic cap (`MAXEXP`) for reproducible counts.
3. Long-square tilting step off the walks; n = 6 c2-3, n = 7 c1+.
4. Row-set comparison n = 15; 4-letter scan n = 16; n = 17 key-coarser lists.
- Watch: a cap is not a verdict (n = 9 zero is weak); parallel timings inflated; equal sizes are not equal sets.
