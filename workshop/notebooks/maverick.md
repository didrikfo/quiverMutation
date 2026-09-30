# Maverick's notebook

## What I now believe (after round 007)
- H-017 (relations > cords among quipus proved in a class outside the quipu theorem) survives depth 4-6 at n = 9 with the ten pairs of the record; per-LNA minimum of rels - cords is depth dependent (4444400: 2 at depth 5, 1 at depth 6); class minimum is what matters.
- The Coxeter polynomial cannot see (cords, relations): c_{n-1} = 1 for every tree algebra with monomial relations; c_{n-2} is not a function of them. Smith form of C+C^T and F-047 profile keep identical survivors, incl. below-diagonal ones.
- Euler signature: pos(C+C^T) <= n-2 iff outside every quipu class for LNAs n = 8..11; forced to fail at n = 13. One relation lowers pos by at most 1, so relations >= 1 or 2, not > cords. (Modest; near F-045/F-048.)
- Round 007 (tested on small cases): the verify-style mutation search (families.verify / linesReachedFrom) DOES return positives: round trip from certified quipu-with-relations members back to their source LNA, 273/273 at n = 7 (L = 4), 84/84 and 126/126 at n = 6; 0/84 at depth L-1. It has resolution exactly its depth. So round 004's depth-4 "reached nothing" on 16 below-diagonal candidates means only "nothing within 4 steps"; forward walks needed depth 6 for 3033030/4444400.

## What I tried
- census of candidates by (cords,rels) at n=9; proved-member walks; Euler/profile filters; search from 16 below-diagonal survivors at depth 4 (round 004); signature by status n=8..11; quipu pos count to n=15; round-trip control (`rounds/007/maverick_control.py`), n=6,7.

## Watch for
- Cords: count cords with m > 0 in `quipuParameters`.
- "Each relation is a rank-2 perturbation" holds only for gldim <= 2.
- A negative at depth d excludes only distance <= d. State the depth against the depth the forward walk needed.
- The n = 7 control run hit the 10 min cap; `maverick_control.py` has no budget flag and prints its summary last.

## Next
- Rerun the 16 below-diagonal candidates at depth 5-6 (size with K = 1; OVERNIGHT proposal if > 10 min).
- Control at n = 9 with a path-5/6 member of 4444400: confirm the flip at L.
- Meet-in-the-middle search from candidate and class side.
- H-017 at n = 10 (262 LNAs, depth 4) and depth 7 at n = 9: overnight. Signature criterion at n = 12, 13.
- Unasked: is the tubular class's corank-2 Euler form a Z-lattice invariant naming its quipu-with-relations members (normal form for H-014)?
