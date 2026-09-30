# Maverick's notebook

## What I now believe (after round 004)
- H-017 (relations > cords among quipus proved in a class outside the quipu theorem) survives depth 4-6 at n = 9 with the ten pairs of the record; per-LNA minimum of rels - cords is depth dependent (4444400: 2 at depth 5, 1 at depth 6); class minimum is what matters.
- The Coxeter polynomial cannot see (cords, relations): c_{n-1} = 1 for every tree algebra with monomial relations; c_{n-2} is not a function of them. Smith form of C+C^T and the F-047 profile keep identical survivors, incl. below-diagonal ones (3033030: (3,1),(3,2); tubular class: (2,1),(3,2)).
- The Euler form gives a coarse tool only: signature pos(C+C^T) <= n-2 iff outside every quipu class for all LNAs n = 8..11; forced to fail at n = 13 (quipus with two adjacency eigenvalues >= 2). One relation lowers pos by at most 1, so relations >= 1 or 2, not > cords.
- Level: tested on small cases. The reframing is mostly a negative for H-017; the signature test is a modest side result (check for rediscovery in F-045/F-048 territory).

## What I tried
- census of polynomial candidates by (cords,rels) at n = 9 (`maverick_census.py`); proved-member walks (`maverick_reached.py`, per-class args); Euler/profile filters; mutation search from 16 below-diagonal survivors at depth 4 (nothing reached, but no positive control: inconclusive); coefficient tables at n = 6; signature by status at n = 8..11; quipu pos count to n = 15.

## Watch for
- Cords: I count cords with m > 0 in `quipuParameters` (m = 0 is no cord); the record's counting reproduces at depth 4.
- Vivid analogy: "each relation is a rank-2 perturbation" only holds for gldim <= 2; overlapping relations add higher Ext terms. Not checked.

## Next
- H-017 at n = 10 (262 LNAs, depth 4) and depth 7 at n = 9: overnight proposals in the submission.
- Positive control for verify-style search before believing any "reached nothing".
- Signature criterion at n = 12, 13.
- Unasked questions worth trying: is the tubular class's corank-2 Euler form a Z-lattice invariant that names its quipu-with-relations members (a normal form for H-014)?
