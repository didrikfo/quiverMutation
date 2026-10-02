# Experimentalist notebook (rewritten each round)

## What I now believe (after round 027)
- Conjecture W (E-107: out-degree 2 reject iff two-term relation p1 b1 = p2 b1 through v, x != 0, x b2 = 0) survives round 027: 0 mismatches on 32 132 out-2 rows (n = 8 c0 61 rejects all W; c2, c3, n = 9 c0: 0 out-2 rejects, W all False). Only class 0 supplies positives, so the converse is tested on one family.
- Parallel-arrow rows (1 135 out-2 with parallel out-arrows, 300 out >= 3): never a reject, W never true. No positive control exists for them.
- Out-degree >= 3: 0 rejects in 3 494 rows over all walks. Out-degree 1 rejects: n = 8 c3 6, n = 9 c0 110 (long-square family, E-103/E-109; unclassified here).
- Earlier: J != 0 <=> tiltingPlus False on all rows (025). T5 long square at 100 % of rejects at n <= 7. n = 9 c0 walk ratio about 2.5 per BFS level, will not close in minutes (026 toolsmith).
- Class sizes: n = 8: 11 classes (2, 8, 18, 20, 26, 52, 80, 128, 128, 130, 266); n = 9: 19 classes, c0 2 algebras.

## What I tried
- 027: `rounds/027/experimentalist_w.py n class budget [maxexp]` (walk + inline W + kerdim, keyed arrows); 4 parallel 540 s runs (about 18 expansions/s each, vs 40 solo): n = 8 c0, c2, c3, n = 9 c0.
- 025: capped 500 s walks; parallel cover fewer algebras than solo. Quote rates, not counts.

## What I would do next
1. Hand-built positive controls for W: parallel b1, and the "cancels" branch (x b2 = 0 only as a sum).
2. n = 8 classes 4-10 (bigger, more varied) for out-degree 2 rejects; solo runs with maxexp.
3. n = 9 c0 overnight with `toolsmith_rejwalk.py` plus W tally (proposal in 027 submission).
4. Classify the out-1 J != 0 rows at c3/n = 9 by long square (E-109 `longSquare`).
- Watch: a cap is not a verdict (zeros at c2, c3, n = 9 are agreement of negatives); parallel timings inflated; seen != expanded.
