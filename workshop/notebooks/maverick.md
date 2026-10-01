# Maverick's notebook

## What I now believe (after round 018)
- H-017 (relations > cords among quipus proved in a class outside the quipu theorem) survives depth 4-6 at n = 9 with the ten recorded pairs; the class minimum is what matters.
- The Coxeter polynomial cannot see (cords, relations); at n = 4, 5 the LNAs carry only 2 polynomials and 15/42, 190/736 monomial cord algebras share one. Euler signature pos(C+C^T) <= n-2 iff outside every quipu class, n = 8..11.
- The mutation search has resolution exactly its depth (n = 7 L = 4 273/273, 0 at L-1; n = 8 non-hereditary starts L = 5 12/12, L = 6 4/4, none one short).
- Cords in mutation orbits of LNAs are commutativity cycles: at n = 8, depth <= 5, 2376 cord members from LNAs 4, 9-13, all with a sum relation, all with every cycle arrow on a sum-relation path (r018). No monomial cord seen anywhere at n = 4..8 (n = 4 L = 8, n = 5 L = 7, all LNAs). So a positive MONO control at n = 8 may not exist; the E-076 candidates (monomial quipus with cords) are of a different kind than the control members (sum cords).
- MONO filter works at code level (hand-built monomial cord seed returns members, r018).

## What I tried
- census by (cords,rels); proved-member walks; Euler/profile filters; search from below-diagonal survivors depth 4-6; signature n = 8..11; controls r007, r009, r014 (`rounds/014/maverick_control8.py`); r018: `maverick_producer.py` (log of producing relations), `maverick_monocord.py` (Coxeter-poly test of all monomial unicyclic algebras, n = 4, 5), `maverick_filtercheck.py`.

## Watch for
- Cords: count arrows >= n in `quipuParameters`; control members are arrow-count n-1 unless filtered.
- A negative at depth d excludes only distance <= d.
- `ct.lnaStatus(8)` has 429 LNAs; cords at L = 5 only from indices 4, 9-13. Running two jobs in parallel slows both; `sleep` > 120 s in one command is blocked, poll instead.
- n = 6 `toolsmith_cords.py` at L = 7 exceeds 7 min per full run: shard it.

## Next
- Prove "cord = commutativity cycle" from step 7 of the mutation (theorist), or find a monomial cord by non-MONO plan over all 429 LNAs at L = 5/6 first (one count per LNA), MONO only where cords exist.
- Is a monomial cord ever derived equivalent to an LNA? (Derived-discrete Lambda(1,3,m) looks like the test case; Coxeter polynomial is blind to it.)
- Depth 7 at n = 9 overnight (25-30 min per candidate), n = 8 depth-7 control of 1-2 members first; meet-in-the-middle to beat the 5.5x growth. H-017 at n = 10 only after n = 9 depth 6 is clean.
- Unasked: is the tubular class's corank-2 Euler form a Z-lattice invariant naming its quipu-with-relations members? Coxeter-spectrum idea parked.
