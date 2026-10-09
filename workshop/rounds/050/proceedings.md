# Round 050 -- proceedings (ordinary)

Decisions from round 049's questions (step 0.5): q1 (overnight depth-7 child ball) -- no, the toolsmith sharded it under 10 min; q2 (PDFs) -- cannot be changed from inside the workshop. Both recorded in STEERING.md.

## toolsmith -- depth-7 child ball for the E-155 misses (T10 i)
Claim: under the J = 0 premise, c1 child 12 and 3 of the 5 misses are joined at total <= 12; c1 14/15 (one key) stay open. Referee (theorist): minor revision (count "4 of 5 misses" wrong, no-key counts missing, E-153 already records no-key nodes, premise in title). Response did items 1-5; item 6 (control on the no-key route) deferred. Decision: **accept, narrowed** (the open control is stated in the entry). Promoted: **E-158**.

## skeptic -- independent Hom(T,T[m]) test of the E-155 premise
Claim: all 40 printed paths (324 edges) pass; the 25 failing J != 0 steps are rejected. Referee (experimentalist): minor revision; every count reproduced; wording only. The chair made the two wording fixes (title "generation assumed"; Cartan vs Cartan^T). Decision: **accept**. Promoted: **E-159**. Referee's point kept: the test is Ladkani's criterion restated, so it adds code independence, not new information.

## maverick -- H-017 Euler signature (T6, dormant since round 023)
Claim: the signature is Cartan-determined, hence a class invariant, blind to (cords, relations); 16 new n = 11 separations. Referee (scholar): minor revision (prior record E-063/E-152; F-048 overlap; shared polynomials). Response did all four items; F-010 Smith-form check not done. Decision: **accept**. Promoted: **E-160**. H-017 status unchanged (OPEN).

## Consequences
- T10: 23 of 25 failing children plus the 25 parents are joined to an LNA by J = 0 paths whose steps passed an independent test (E-155 paths: E-159; the 3 new paths: not yet replayed with the independent test). Open: c1 14/15 and the new paths. No claim yet that the guard's failures stay in the class.
- `fingerprint.canonicalKey` docstring ("no search has yet produced" None) is wrong (E-158). Not edited: a library docstring change waits for a toolsmith round with the T10 (iv) reword.
- H-017: the Euler signature cannot be used as a test of its relations-vs-cords claim (E-160); thread T6 closed as to the signature, H-017 stays OPEN.

## Questions for the steering committee
1. **Overnight: c1 14/15 at horizon 13** (about 2 h on 4 cores, or a depth-6 target ball ~20 min first). Recommend: no; toolsmith first tries a target ball of depth 6 (shardable, under 10 min per shard) and a larger `canonicalKey` cap. (Unanswered, the next chair keeps that.)
2. **PDFs of arXiv:1009.3370 and 2509.12983**, if you can supply them. Recommend: yes.

## Decisions taken for the steering committee
- Round 049 q1: no overnight; sharded instead (done, E-158). Round 049 q2: unchanged, literature parked.
