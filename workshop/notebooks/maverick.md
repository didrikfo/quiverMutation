# Maverick's notebook

## What I now believe (after round 009)
- H-017 (relations > cords among quipus proved in a class outside the quipu theorem) survives depth 4-6 at n = 9 with the ten recorded pairs; per-LNA minimum depends on depth (4444400: 2 at 5, 1 at 6); the class minimum is what matters.
- The Coxeter polynomial cannot see (cords, relations) (c_{n-1} = 1 for monomial trees). Euler signature pos(C+C^T) <= n-2 iff outside every quipu class, n = 8..11; one relation lowers pos by at most 1, so relations >= 1 or 2, not > cords.
- The mutation search (families.verify / linesReachedFrom) has resolution exactly its depth: round trips 273/273 at n = 7 (L = 4), 84/84 at n = 6 (L=3); round 009: L = 5 passes 42/42 at n = 6, 8/8 at n = 7, 0/50 at L-1; recorded paths are shortest (reachedQuipuAlgebras is an exhaustive DFS keeping the min).
- Depth-5 search from the four K = 1 below-diagonal candidates at n = 9 reached nothing (new, round 009); depth 6 for candidate 1 also nothing. So candidate negatives now exclude members within 5 (one: 6) steps; forward walks needed 6.
- Cost per depth level at n = 9 is about 5.5-5.7x: depth 4 ~14 s, 5 ~80 s, 6 ~5-7 min, 7 ~27 min per candidate.

## What I tried
- census by (cords,rels); proved-member walks; Euler/profile filters; search from below-diagonal survivors at depth 4 (r004), 5 and 6 (r009); signature n=8..11; control rounds 007 (L=3,4), 009 (L=5, shortest check; `rounds/009/maverick_control5.py`, has --plan).

## Watch for
- Cords: count cords with m > 0 in `quipuParameters`.
- A negative at depth d excludes only distance <= d; say it next to the depth the forward walk needed (6).
- `ct.lnaStatus(6)` has 42 LNAs, not 84 (E-069's 84 counts members).
- n = 7 full control (132 LNAs, ~14 s each) exceeds 10 min; shard it.
- `maverick_verify.py` has no candidate index: K = 1 gives 4 cells, K larger gives the 16.

## Next
- Depth 6 on all 16 candidates: ~1.5 h, shards of 2 per command (needs a candidate-index arg). Depth 7: ~7 h, overnight, already in Menu 4 for the LNA side.
- Meet-in-the-middle from candidate and class side to beat the 5.7x growth.
- H-017 at n = 10 (262 LNAs) only once depth 6 at n = 9 is clean. Signature criterion at n = 12, 13.
- Unasked: is the tubular class's corank-2 Euler form a Z-lattice invariant naming its quipu-with-relations members? Coxeter-spectrum idea parked.
