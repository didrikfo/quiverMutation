# Theorist position (round 012)

## Most promising question
**Can move sequences within single-relation orbits P, Q characterize parity-class words through their internal structure alone—without fitting to the lists A, B?**

Mechanism over pattern: E-077 shows A/B alternate between P and Q, but the rule is fitted to data. The deeper question is whether the mutation rule's move graph, restricted to single-relation rows of keys `35`/`36`, has an invariant (reachability, centrality, orbit-pair alignment) that forces this alternation. If yes, the pattern is not accidental and survives n > 20.

## Weakest claim
That no algebraic structure (GF(2) functional, SNF, integer mod 2/4) separates P from Q (notebook, round 011). This was checked against three specific methods; absence-of-evidence is not evidence-of-absence. Worse: I assumed those methods were the only natural ones, without trying the move-graph structure itself.

## What I need
1. **From experimentalist/toolsmith:** the full BFS move sequences connecting `35@0 ↔ 35@2` in P at n = 12 (have 178 rows, need 1–2 edge moves); same for `4@0 ↔ 3@1` across P and Q; show that `5@0` has no such move into `3@2` (what parity obstruction blocks it?).
2. **From experimentalist:** structure of `5046`, `5056` orbits at n = 17 (or `--plan` one early in round 013 if too large); do they live in a third single-relation key class, or do they split like `4056` (orbit vs. key)?
3. **From scholar:** do the 7 cores of E-059 carry single-relation structure in any key? (quick check: are their sorted-key tables sparse or dense?)

**Why this unblocks work:** once we know the move graph's role, the 7 cores either fit or don't, and we can rule out spectral methods altogether.
