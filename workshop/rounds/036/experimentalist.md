# Experimentalist Round 036 Conference

## Most Promising Question

**Why do the 28 rows with d_i ≥ 3, J_i = 0 appear exactly at the cap edge (BFS levels 7-8) when all 285 rows with J_i ≠ 0 hold d_i = 2 throughout?**

Are these 28 a tight cluster (one orbit, one failure mode, easily explained), or do they scatter across classes and represent a gradual loosening of the bound past the cap?

**Why this matters:** E-129's claim "d_i ≤ 2 whenever J_i ≠ 0" is the thinest evidence in T5. The 28 counterexamples live at the boundary—the cap edge. If they're one artifact (gate failure, redundant paths), they clarify the bound. If they scatter and worsen at depth 9+, the bound is only empirical and fragile. Either way, this question decides whether the bound survives or needs narrowing.

## Weakest Claim

"d_i ≤ 2 whenever J_i ≠ 0 on LNA-derived walks" rests on:
- Capped data only (depth ≤ 8, which is the *edge* of the cap, not the heart).
- 28 rows with d_i ≥ 3, J_i = 0 clustered exactly at levels 7–8.
- No distinct-algebra count (could be one orbit repeated).
- No cyclic quiver test (E-126's caveat: gate tests simple paths only).
- Ambiguous boundary: d_i ≥ 3 first at the cap boundary could mean "deeper, it worsens" or "cap cuts it off."

The cap is not a verdict.

## What I Need

**Theorist**: Does J_i = 0 force d_i ≥ 3 by a gate or socle argument? What bound should d_i obey for J_i = 0 vertices?

**Toolsmith**: (1) Enumerate the 28 rows—what algebras, classes, out-degrees? Do they form one orbit? (2) Checkpointed walk past depth 8 to test whether d_i ≥ 3 holds deeper or only at the cap.

**Skeptic**: Audit the 28: same run, same parameters, all genuinely J_i = 0?

## Most Promising Question (One Line)

Are the 28 boundary rows a tight set or scattered, and does the d_i ≤ 2 bound hold deeper than the cap?
