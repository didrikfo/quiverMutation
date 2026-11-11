# Experimentalist, round 124 (phase 1 conference, no new work)

## Most promising question for the next phase
Which K = 0 rows transit between K >= 1 sub-orbits carrying different image keys (H-026 / S-1)? The
round-123 review showed that constancy on K >= 1 sub-orbits is automatic (moves induce image moves,
E-141, E-146) and the image key can change only across rows with a relation at vertex 1. So the class-level
statement is a statement about K = 0 transit. It is cheap: the data (orbits, paths, scripts) exist at
n = 13 and n = 15. It would also explain the `rules` contrast (n = 15 h = 3 constant, n = 13 h = 2 mixed).
## One promising question outside the thread with most of this phase
T7 / H-010 (dormant since 119, E-313 had no power): fix the engine ValueError on parallel arrows
(liftPath at k = 5, 6) and run start-gap-2 pairs with a planted lowering as power control. Without it, the
retirement route for H-010 rests on a test that cannot say YES. Needs a toolsmith patch and an overnight slot.
## Weakest claim the workshop relies on
That every row-to-row step in the move tables (`rules`, `edgeMoves`, `doubleMutation`) is a derived
equivalence. E-141, E-146, my E-307 extension and the round-123 paths all assume it; maverick checked the
paths against the tables, not against the mutation engine. E-318 is a reminder that a guarded step can look
admissible and not be one. A second weak one: "ten odd-n fits" for K0 = min(h, n-4-h) are fits, with
shuffled-label controls only at n = 13, 15.
## Result closest to a finding
E-318: a key-kept J != 0 right mutation at n = 11 changes SNF(S) tail (2,4) to (2,2), so the guarded walk
leaves its derived class (H-015 refuted in its step reading). One check still needed: an independent
rebuild by a non-author of the step and of the SNF(S) of both endpoints (and that both are acyclic, det C = 1),
since it is one hit in 20 roots and the whole T10 conclusion rests on it. E-321 (no library merge uses a
J != 0 step) is the natural companion, already a replay.
## What I need
- theorist: the right class-level H-026 statement (K = 0 transit, or image key set per orbit) before I run more.
- toolsmith: mutation-engine replay of the 26-move n = 15 and 13-move n = 13 h = 2 paths; the ValueError
  patch for T7; a faster Coxeter key for n = 17 (2.67 M rows at n = 15 already takes 165 s).
- skeptic: independent rebuild of E-318 (see above).
Owed: my revision of rounds/123/experimentalist.md per its review.
