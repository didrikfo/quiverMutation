# Toolsmith round 020

## Most Promising Question

**Can we predict from LNA structure which algebras yield cord members, eliminating wasted depth-6 walks on certain negatives?**

E-089 and E-076 together suggest cord membership is determinable. E-089 found no monomial cord at n=4,5,8 across six LNAs, but E-076's candidates are monomial with cords. If LNA shape or dimension correlates with cord existence, we can shard the 429 remaining n=8 LNAs by predicted sign before walking—cutting the walk cost by orders of magnitude.

## Weakest Claim

E-087: no monomial cord at n=8. Rests on 2 found members across only 6 walked LNAs (L=6, depth ≤ 5), yet E-076 proves monomial cords exist in the catalogue. The heuristic "a cord needs sum relation" is observed, not derived; E-089 walked at shallow depth. Negatives at n=4,5 are shallow (depth 6); the search is incomplete.

## What I Need

- **Theorist:** algebraic predicate—kernel dimension, radical structure, Euler form—separating sum-relation cords from monomial ones in the output of `reachedQuipuAlgebras`.
- **Experimentalist:** shard count for L=5 MONO --plan over all 429 n=8 LNAs. Thirty minutes per 50 LNAs, tells us whether depth 5 negatives dominate.
