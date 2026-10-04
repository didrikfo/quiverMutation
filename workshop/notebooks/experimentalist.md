# Experimentalist notebook (rewritten each round)

## What I now believe (after round 033)
- A both-die square (p1 b1 = p2 b1, p1 b2 and p2 b2 both monomially killed) is a legal gate-admitted J != 0 vertex at n = 6. Whether a walk can contain it is a question about the Coxeter key: over 4 core shapes + pendants, 0 of 42 such algebras at n = 6 have an LNA key; 48 of them at n = 7 (cores with two length-2 sides, or (2,3)); 1 408 at n = 8. Squares with an arrow side (E-120's n = 6, 7 shapes) never get an LNA key. So below n = 7 the absence is the algebras (key-level, necessary test); at n = 7 it is unresolved: capped BFS (about 37 000 expanded) found 0 of 44 candidates.
- dim J_i <= 1 on every walk J != 0 row at n = 6, 7, 8 (c0, c1), always inside dim e_iAe_v = 2; 1 to 3 vertices i per row. One out-degree 3 row at n = 8 c0 with parallel arrows (conflicts with E-111's "0 out >= 3 rejects"?, unverified).
- Earlier (030): dim e_iAe_v at n = 8 walks grows with depth (8 at depth 8); parallel arrows give dim 2; W has 0 mismatches on 32 132 out-2 rows (027); J != 0 <=> tiltingPlus False (025).
- Class sizes n = 8: 11 classes (2, 8, 18, 20, 26, 52, 80, 128, 128, 130, 266); n = 9: 19 classes.

## What I tried
- 033: `rounds/033/experimentalist_bothdie.py` (hand, enum, reach), `experimentalist_dimji.py`; 4 parallel jobs slow each other about 2x.
- 030: BFS depth tables of max dim; 027: W test walks; 025: capped walks.

## What I would do next
1. Decide the n = 7 case: close the BFS of classes idx 1, 3 (overnight) or reverse-mutate the 44 candidates toward an LNA.
2. Check the out-degree 3 parallel row against E-111; reconcile 136 vs 61 out-2 rows at n = 8 c0.
3. Enumerate both-die cores more widely (two-term kills, two-arrow pendants, scalar != 1) before trusting "arrow-side squares never get an LNA key".
4. n = 9 c0 prefix dim J_i; positive control for W (parallel b1, "cancels").
- Watch: a cap is not a verdict; an LNA key is a necessary test only; parallel runs inflate timings; seen != expanded.
