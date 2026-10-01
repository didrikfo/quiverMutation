# Maverick's notebook

## What I now believe (after round 014)
- H-017 (relations > cords among quipus proved in a class outside the quipu theorem) survives depth 4-6 at n = 9 with the ten recorded pairs; the class minimum is what matters.
- The Coxeter polynomial cannot see (cords, relations). Euler signature pos(C+C^T) <= n-2 iff outside every quipu class, n = 8..11; one relation lowers pos by at most 1, so relations >= 1 or 2, not > cords.
- The mutation search has resolution exactly its depth (n = 7: 273/273 at L = 4, 0 at L-1; n = 6 L = 5 42/42; round 013 n = 7 L = 6 16/16).
- Round 014: n = 8 control from non-hereditary starts: L = 5 12/12 found, 0/12 at depth 4; L = 6 4/4 found (3 LNAs; one 4-relation member), 0/1 at depth 5. Nodes at depth 6: 6.9e3-3.85e4 (n = 9 negatives 5e4-6e4). So the machinery handles relation-bearing starts; but all control members were 7 arrows on 8 vertices (zero cords), so no control with cords > 0 yet.
- Cost: n = 8 depth-6 search < 100 s per member; the walk that finds a member costs 90-310 s per LNA (depth 6), the real bottleneck. Per depth level at n = 9 about 5.5x.

## What I tried
- census by (cords,rels); proved-member walks; Euler/profile filters; search from below-diagonal survivors at depth 4, 5, 6; signature n = 8..11; controls r007 (L=3,4), r009 (L=5), r014 (n=8, nonhereditary; `rounds/014/maverick_control8.py`, has --plan, SHORT=1, HIGH=1).

## Watch for
- Cords: count cords with m > 0 in `quipuParameters`; control members are arrow-count n-1 unless filtered.
- A negative at depth d excludes only distance <= d; the forward walk needed 6.
- `ct.lnaStatus(6)` has 42 LNAs, `(8)` 429. n = 7 full control exceeds 10 min; shard it.
- Mis-typed positional args waste a run (`_L6_plan.txt` is L = 1 by mistake).
- Parallel jobs slow walks 2x; a 10-minute command needs one walk plus one search.

## Next
- Control with cords > 0 and relations >= 1 at n = 8 (filter arrows >= n in member list).
- Depth 7 at n = 9 overnight (about 25-30 min per candidate); an n = 8 depth-7 control of 1-2 members first.
- Meet-in-the-middle from candidate and class side to beat the 5.5x growth.
- H-017 at n = 10 only after depth 6 at n = 9 is clean. Signature criterion n = 12, 13.
- Unasked: is the tubular class's corank-2 Euler form a Z-lattice invariant naming its quipu-with-relations members? Coxeter-spectrum idea parked.
