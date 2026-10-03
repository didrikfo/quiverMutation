# Experimentalist notebook (rewritten each round)

## What I now believe (after round 030)
- dim e_iAe_v (all ordered pairs) at n = 8 c0/c1 walks: <= 1 for depth <= 3, first 2 at depth 4, 3 and 4 at depth 6, 8 at depth 8 (c0). Not bounded by 1 or 2 (round 030, `rounds/030/experimentalist_dimdepth.py`). Parallel arrows alone give dim 2 (first c0 example has a doubled relation), so dim >= 2 is not evidence of a circuit.
- Over out-degree-2 mutation vertices: max 1 to depth 4, 2 from depth 5, 3 from depth 7/8, 4 once (c1, depth 7).
- Conjecture W (E-107) survives 027: 0 mismatches on 32 132 out-2 rows; positives only at n = 8 c0 (61 rejects). No positive control for parallel arrows. Out-degree >= 3: 0 rejects in 3 494 rows. Out-1 rejects: n = 8 c3 6, n = 9 c0 110 (long square, E-103/E-109).
- J != 0 <=> tiltingPlus False on all rows (025). n = 9 c0 walk ratio about 2.5 per BFS level; does not close in minutes.
- Class sizes n = 8: 11 classes (2, 8, 18, 20, 26, 52, 80, 128, 128, 130, 266); n = 9: 19 classes, c0 2 algebras.

## What I tried
- 030: BFS depth tables of max dim, 480 s capped walks c0 (depth 9, 7 679 measured) and c1 (depth 7, 12 430), run in parallel (slower per run).
- 027: `rounds/027/experimentalist_w.py` (walk + W + kerdim); 4 parallel 540 s runs.
- 025: capped 500 s walks; quote rates, not counts.

## What I would do next
1. Pairs with two nonzero non-parallel classes only; Gamma_i shapes at the dim >= 3 out-2 rows.
2. Verify per-level completeness (levels <= 6 whole); n = 9 c0 to depth 6.
3. Positive controls for W: parallel b1, "cancels" branch (hand-built).
4. n = 8 classes 4-10 for out-degree 2 rejects; classify out-1 J != 0 rows by long square.
- Watch: a cap is not a verdict; last BFS level is partial; parallel runs inflate timings; seen != expanded.
