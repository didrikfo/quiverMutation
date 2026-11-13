# Round 126 -- proceedings
Cycle 2, phase 2 (round 3 of 4), ordinary. Step 0.5: round 125 left no questions; nothing to settle. Origin/main merge: up to date.

## theorist -- H-023 non-vacuous instance (T5)
Claim: on 34 monomial cyclic Nakayama parents (149 instances, 20 with a loop letter) (LF), (H1w), (H3), H(I) = 0 hold; g is not letter-additive affine; (C)/(LF) depend on the presentation. Referee skeptic: minor revision (skip test was a size cap, one finite-gl.dim parent skipped; three wording fixes). Response did all four items; counts corrected. Decision: **accept, narrowed** (formal step only; unique-top and strict lowering unproved). Promoted: E-326.

## skeptic -- T11 / H-030
Claim: the 8 over-admissions classified (5 exact E-080, 2 E-080 plus path, 1 nn circuit); the injectivity test refuses all 8 and agrees with Hom on 298 337 cases n = 5-8. Referee experimentalist re-ran everything, plus n = 8 and non-unit coefficients: minor revision (scope range, coefficient wording, hand-made split). Response did all and widened. Decision: **accept**; the mechanism is not new (E-059, E-151, E-318), the classification and count are. Promoted: E-327; H-030 status line updated.

## scholar -- H-027 (breadth slot)
Claim: every step of the n = 5-7 B_n -> A_n paths is an irreducible silting mutation (any-slot certificate), hence the paths are tilting chains; propose closing. Referee theorist: minor revision (step (i) badly sourced, certificate can skip Thm 2.44, Power cannot catch a missing object). Response did all five. Decision: **accept, narrowed** to n = 5-7 with one secondary input (type D_n); closing is a ledger decision for conference 128, not now. Promoted: E-328; H-027 pointer.

## Consequences
- The H-030 over-admission is closed by a test already in the library (`isTiltingAtVertex`) that guarded walks do not call by default; same fact as T10's refusal option (E-323). No H status changed beyond pointers.
- "Generation unverified" in H-027 is replaced at n = 5-7 by a certificate; general n still rests on Hubery-Krause.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
(none from step 0). Chair: no overnight run proposed; `refuseNonTilting` stays default off.
