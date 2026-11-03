# Round 111 -- proceedings
Cycle 1, phase 3 (research), ordinary. Merge of origin/main: already up to date. Step 0.5: round 110 left no questions; nothing to decide.

## toolsmith -- T10: n = 7 class 1, the key-kept J != 0 steps (referee scholar: minor revision, answered)
Claim (after response): all 89 key-kept J != 0, D = 0 steps of the 45 000-expansion class-1 walk (E-201's 13 included) are Z-congruent parent to child by an exhaustive search, and so are E-147's 4 pickled class-2 steps (by exhibit); this is the expected same-key result (E-242) and the test has no NO case at n = 7. Referee: YES reproduced 89 of 89; asked for E-242 framing, reconciliation with E-147's "incongruent", step identifiers. All done. **Accept, narrowed**: the walk (25 min) was not re-run by the referee; class 2 is 4 of 9 steps. Promoted: E-295. Does not unblock the R-candidate; it shows the Cartan route cannot.

## theorist -- T5: loop model M2 (referee skeptic: minor revision, answered)
Claim: the cancellation of Prop. 7.6 for the formal dg algebra k<a,b,f> (free degree-1 loop) is proved by a g-filtration reducing to Lemma A; a degree-0 loop with weight 3 also cancels to weight 12. Referee: proof sound and reproduced, new; softer wording for the degree-0 part (Lemma A tested for degrees 0..2 only), quote E-284. Response ran Lemma A with degree -1 arrows (0 mismatches, weight <= 6) and narrowed the title. I re-ran two scripts: same output. **Accept, narrowed**: M2 is a formal model, no Hom/End(T) computation, M4 and M5 untouched. Promoted: E-294. H-023 stays OPEN.

## maverick -- S-1: K0 = 6 at n = 17 (referee theorist: minor revision, answered)
Claim: at key level, head ends, one lone-3 key class per n, the K-threshold law holds at n = 17 (K0 = 6) and n = 19 (K0 = 7); tail ends by tested mirror equivariance. Referee reproduced n = 15 and n = 17 K1 = 6; found errors in the "in key class" column and an undefined K0 convention. Response fixed the n = 15 count (527), stated the convention, ran the mirror test (0 failures) and narrowed the scope. **Accept, narrowed**: n = 19 and the n = 17 K = 5 cell were not re-run by anyone. Promoted: E-296; H-026 status line extended (stays OPEN).

## Consequences
- E-284's remark "degree-0 loops are not testable" is too strong at w(g) = 3 (E-294); no status changes, but the T5 degree-0 hypothesis claim now has a computation behind it, in one model.
- H-015 / R-candidate "a guarded walk stays in one class": the Cartan-level test cannot separate key-kept steps at n = 7 (E-295); the check that blocks the R-entry is now a Hom-verified tilting chain for the 13 steps of E-201 (parents with parallel arrows need the Hom test; cf. E-212 record 6), or an invariant finer than the Cartan matrix. Thread T10 updated.
- F-candidate H-024: no change (its blocking check is the lam < 1.5 NO side).

## Questions for the steering committee
None.

## Decisions taken for the steering committee
Step 0.5: nothing to decide. No overnight run approved. The n = 19 re-run of E-296 (`maverick_k17.py run 19 7`, over 10 minutes in parallel, about 9 minutes alone) is left to an ordinary round.
