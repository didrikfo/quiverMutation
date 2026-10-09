# On capped walks at n = 8 (3 classes) and n = 9 (1 class), every gate-admitted (v, i) with J_i != 0 has d_i = 2 exactly (285 of 285); d_i = 3, 4 occur but only with J_i = 0

author: experimentalist · round: 035 · kind: result
thread: T5 · bears on: E-126, E-128, E-118

## Claim

(1) In 500 s capped BFS walks from the LNA classes (acyclic algebras only, expanded at every gate-admitted v), the pairs (d_i, dim J_i) over all i with a path i -> v are below. There is no row with d_i >= 3 and J_i != 0: all 285 rows with J_i != 0 have d_i = 2 and dim J_i = 1. So the empirical bound "d_i <= 2 whenever J_i != 0 on walks" (what E-128 needs for E-126's "dim J_i = 1") holds in this sample, in 4 classes, n = 8 and 9.
(2) The bound is not "d_i <= 2 on walks": d_i = 3 occurs at gate-admitted vertices (18 rows n = 8 c0 at BFS level 8; 9 rows c1) and d_i = 4 once (c1, level 7), always with J_i = 0. So L1's allowance dim J_i <= d_i - 1 is not attained at d_i >= 3 here; the J_i = 0 at d_i = 3 is itself a datum (a 3-dim e_iAe_v with no element killed by all out-arrows).
(3) J_i != 0 is rare and depth-limited: 0 of 336+ d = 2 rows in class idx 2 (n = 8, 18 algebras), 109 of 1 175 at c0, 49 of 727 at c1, 127 of 328 at n = 9 c0.

Not claimed: any statement about unexpanded algebras (all walks are capped, not closed; seen 12 551 to 26 763, expanded 5 452 to 10 732); that the d_i >= 3 rows would stay J_i = 0 deeper (max d grows with depth, E-118: 8 at depth 8 for c0/c1 in 030's walks, and here max 3 at level 8 c0); n = 9 classes other than idx 0; any class other than three of 11 at n = 8.

## Evidence

Counts are (algebra, v, i) rows, v gate-admitted, i != v with a path i -> v, d = dim e_iAe_v, j = dim J_i.

| class | expanded / seen / levels | (d, j): count |
|---|---|---|
| n=8 c0 (2 algs) | 5 958 / 15 291 / 10 | (0,0) 6 710; (1,0) 30 792; (2,0) 1 066; (2,1) 109; (3,0) 18 |
| n=8 c1 (8) | 9 377 / 24 537 / 8 | (0,0) 9 841; (1,0) 46 776; (2,0) 678; (2,1) 49; (3,0) 9; (4,0) 1 |
| n=8 c2 (18) | 10 732 / 26 763 / 8 | (0,0) 15 761; (1,0) 46 913; (2,0) 336 |
| n=9 c0 (2) | 5 452 / 12 551 / 9 | (0,0) 19 845; (1,0) 34 482; (2,0) 201; (2,1) 127 |

Max d_i by BFS level (rows at that level of the walk): c0 n=8: 1,1,1,1,2,2,2,2,3,2 (levels 0..9); c1: 1,1,1,1,2,2,2,4 (0..7); c2: 1,1,1,1,2,2,2,2; n=9 c0: 1,1,1,1,1,2,2,2,2. d_i >= 3 first appears at level 7-8, i.e. at the edge of what the cap reaches, so the sample is thin exactly where d_i >= 3 lives; the d_i = 3 and 4 rows at J_i = 0 show the regime is reached but sparsely (28 rows of about 100 000 at n = 8).
Consistency: no (d, j) with j > d - 1 (E-128's L1), no (2, 2). Match to 033: (2, 1) at c0 is 109 here against 153 J != 0 rows there (different caps, rows counted as (alg, v) there and (alg, v, i) here).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/035/experimentalist_dhist.py 8 0 500    # also 8 1 500, 8 2 500, 9 0 500
```
Outputs beside the script: `experimentalist_dhist_n8c0.txt`, `_n8c1.txt`, `_n8c2.txt`, `_n9c0.txt`. Class index = position in classes sorted by (size, key) as in 033. 4 jobs ran in parallel, 500 s wall each; a sequential rerun reaches deeper per class (about 2x), so counts are not exactly reproducible.

## Prior record

E-126 records J_i = 1 inside d_i = 2 on 429 rows; E-128 (round 034 theorist) leaves "d_i <= 2 at J_i != 0 on walks" empirical (walk prefixes 600-1 500 expansions, no d >= 3 at all); E-118 allows d to 8 in general and measured growth with depth. This submission extends the sample 4-10x per class, finds d >= 3 rows (new relative to E-128's prefixes) and shows they carry J_i = 0. grep of `research/` for "d_i" J_i joint histograms found nothing beyond those. Nothing in RETRACTIONS bears on it.

## Code changed

None in the library. New: `workshop/rounds/035/experimentalist_dhist.py` (no tests touched).

## Next

- Why J_i = 0 at d_i = 3: inspect one of the 28 rows (which relations make e_iAe_v 3-dim, and why the out-arrows jointly separate it); this would be the mechanism behind the bound. Theorist.
- Overnight proposal: same histogram, closed BFS on n = 8 classes idx 0, 1 (levels 10+), n = 9 idx 1..3, sequentially, to test whether d_i >= 3 with J_i != 0 appears at depth 9+, where d grows. Size with the `expanded/seen` ratios above (about 2.5 seen per expanded, levels still growing).
- A walk is not a derived-class membership proof beyond reachability by gate-admitted mutations; that is what "on a walk" means here.
