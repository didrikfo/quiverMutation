# Maverick's notebook

## What I now believe (after round 021)
- H-017 (relations > cords among quipus proved in a class outside the quipu theorem) survives depth 4-6 at n = 9 with the ten recorded pairs; the class minimum is what matters.
- The Coxeter polynomial cannot see (cords, relations); Euler signature pos(C+C^T) <= n-2 iff outside every quipu class, n = 8..11.
- The mutation search has resolution exactly its depth (n = 7 L = 4 273/273; n = 8 L = 5 12/12, L = 6 4/4, none one short).
- Cords in mutation orbits of LNAs are commutativity cycles (r018): n = 8, depth <= 5, 2376 members, all with a sum relation. No monomial cord at n = 4..8.
- r021: an n = 8 LNA has cord members within 3 steps iff some relation has >= 3 arrows (365/429, 0 mismatches; the 64 rad^2-zero LNAs have none at L = 3 or 5). MONO at L = 3: 0 of 429. Unproved; heuristic is that a >= 3-arrow monomial turns into a commutativity square.
- Sizing: non-MONO L = 5 over the other 365 LNAs is about 4.5 CPU-hours (30-60 s each), less valuable now that C says who has cords.

## What I tried
- census by (cords,rels); proved-member walks; Euler/profile filters; controls r007-r014; r018 producer/monocord/filtercheck; r021 `maverick_predict.py` (cord count per LNA, IDX= and MONO= env), `maverick_criteria.py`.

## Watch for
- Cords: arrows >= n in the raw visitor; `reachedQuipuAlgebras` keeps monomial quipu trees only.
- A negative at depth d excludes only distance <= d. The 64 negatives are at L = 5.
- Parallel jobs slow each other; `sleep` > 120 s in one command is blocked, poll.
- Zero-cord LNAs walk at ~14 s for L = 5, cord LNAs 30-60 s.

## Next
- Prove criterion C (theorist), test at n = 9 (L = 3, cheap), n = 6, 7 for sanity.
- A monomial cord, if any, is not near LNAs: try MONO at L = 5 on the 365 in shards (overnight) or look at derived-discrete Lambda(1,3,m).
- Depth 7 at n = 9 overnight remains for H-017; n = 10 only after n = 9 depth 6 is clean.
- Unasked: is the tubular class's corank-2 Euler form a Z-lattice invariant naming its quipu-with-relations members?
