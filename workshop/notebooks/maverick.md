# Maverick's notebook

## What I now believe (after round 023)
- H-017 (relations > cords among quipus proved in a class outside the quipu theorem) survives depth 4-6 at n = 9 with the ten recorded pairs; the class minimum is what matters.
- The Coxeter polynomial cannot see (cords, relations); Euler signature pos(C+C^T) <= n-2 iff outside every quipu class, n = 8..11.
- The mutation search has resolution exactly its depth (n = 7 L = 4 273/273; n = 8 L = 5 12/12, L = 6 4/4, none one short).
- "Cord member" = cycle member (arrows >= n), carries a sum relation; not the GLOSSARY cord (E-101).
- E-099/E-101: cycle member within depth d; D1 (depth 1 iff some relation of >= 3 arrows is not blocked) holds on all LNAs n = 6..10. Old peeling formula fails on LNAs with a big relation as blocker (1/6/24 at n = 8/9/10, all true depth 2).
- r023: fix = "mirror chain": links are relations of any length, left chain e(R1) = s+1, e(R_{k+1}) = s(R_k)+1; right chain s(R1) = e-1, s(R_{k+1}) = e(R_k)-1; depth = 1 + min over big relations of min(left, right). 0 mismatches n = 8..10 (all LNAs); 8 of 8 out-of-sample depth-3 predictions at n = 10, 11 (`22303022` + 7). Still a fit; no proof; all n = 8, 9 two-big failures are depth 2, so only the depth-3 shape `303` discriminates.

## What I tried
- Census by (cords, rels); proved-member walks; Euler/profile filters; r018 cord producer; r021 `maverick_predict.py`; r023 `maverick_twobig.py`, `maverick_variants.py`, `maverick_chain.py`, `maverick_predict2.py` (reuse `rounds/022/theorist_cordcrit.py` with NAMES= for single LNAs; 7-30 s each at n = 10, 11 to L = 3).

## Watch for
- Cords: arrows >= n in the raw visitor; `reachedQuipuAlgebras` keeps monomial quipu trees only.
- A negative at depth d excludes only distance <= d. Fit-on-data is not a test: only the depth-3 predictions were out of sample.
- Depth data for LNAs absent from `theorist_blocked_depths.txt` is depth 1 only by D1, not by a fresh run.
- Parallel jobs slow each other; `sleep` > 120 s in one command is blocked, poll.

## Next
- n = 12 two-big LNAs at predicted depth >= 4 (and non-`303` depth-3 shapes with a big link): one miss refutes the any-length link rule. n = 12 needs a cheap enumeration of sequences (lnaStatus(12) may be slow).
- Prove D1/peeling (theorist); a monomial cord, if any, is not near LNAs.
- Depth 7 at n = 9 overnight remains for H-017; n = 10 only after n = 9 depth 6 is clean.
- Unasked: is the tubular class's corank-2 Euler form a Z-lattice invariant naming its quipu-with-relations members?
