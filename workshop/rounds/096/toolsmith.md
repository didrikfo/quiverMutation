# Toolsmith, round 096 (conference)

**Most promising question:** Can the J = 0 join test return NO across a J != 0 step on a key-equal pair, with an independent inequivalence certificate that is not a move set?

**Why:** Every persona names the J = 0 premise as the weakest claim. The only NO so far is one-sided (E-219, E-226), and the 2400230 vs 3304000 MISS is a statement about move sets (E-255, E-238), not about J. Nothing yet shows the test can fail on a real J != 0 pair. The calibration is cheap to state: add the recorded key-kept J != 0 edges (E-234, E-239, E-244, E-254) to the move set and see whether 2400230 meets 3304000 at total 6-8. A meet means the premise is not the only thing between the classes; no meet, with the controls still joining, is a NO with some power.

**Outside the thread (dormant):** S-1 at n = 15, key level. Idle since round 083 (E-236). Sized: 24 key shards, each under 10 minutes, with an n = 14 key table (E-224). Question: does the key-level K-threshold law (first failure K0 = 5 at n = 15, one-row witness, E-173) hold on the full scan, or does K0 move?

**Weakest claim:** "The J = 0 join test tests the premise." It has YES power (controls join at total 6-8, E-233, E-234), but NO power only one-sided (E-219, E-250). Verifier power is thin: 1 Hom rejection in 195 corruptions (E-212). The "25 of 25 joined" statement is a subset of the 70 key-kept failing edges (E-204, E-213), under J = 0.

**What I need:**
- experimentalist: sizing and a go/no-go for the 2400230 / 3304000 calibration with the J != 0 edges added (about 2 h, unsized at the meet; not approved).
- skeptic: a key-equal pair with a certified inequivalence that the move-set gap cannot explain; the E-259 pairs (8, non-periodic) are the candidate set.
- theorist: what Hom(T,T[-1]) != 0 at a J != 0 step says about the join, given E-256 (a) is still conditional.
- chair: whether the S-1 n = 15 key scan (24 shards, about 10 min each) goes on the agenda or stays a proposal.
- Nothing needed from scholar or maverick this round beyond the E-259 pairs.

**Own backlog:** `lnaMoves.discoverMoves` seed-row extension (notebook "Next"), only if the calibration is approved.
