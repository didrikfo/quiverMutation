# Round 110 -- proceedings
Cycle 1, phase 3 (research), ordinary. Merge of origin/main: already up to date. Step 0.5: round 109 left no questions; nothing to decide.

## skeptic -- H-024: a pair inequivalent by construction (referee toolsmith: minor revision, answered)
Claim: E-269's Z-congruence search says NO on 19 of 19 LNA pairs made inequivalent by q_S (n = 9..12, lam 1.661..2.46) and YES on 19 of 19 same-quipu controls; each NO coincides with the q_S certificate; no two-quipu LNA pair has lam < 1.5. Referee rebuilt q_S independently (differs on all 19; equal on the group of 140) and re-ran n = 9, 11. Required items (14/5 split of strong/weak NOs with slack-5 counts, nsol table, "never a pure lattice-search NO", cite E-265/E-291) all done. **Accept, narrowed**: power of the NO side shown at lam >= 1.661 on q_S-separated pairs; the lam < 1.5 half stays untested; three n = 11 NOs rest on 1-2 vectors. Promoted: E-292; H-024 status line extended (stays SUPPORTED). Not an F-entry: the NO is never a pure lattice-search NO.

## experimentalist -- T10: a J != 0 key-refusal control (referee skeptic: minor revision, answered)
Claim (after rescoping): at n = 6 key class 0 all 1 150 gate-admitted J != 0 steps of a capped 480 s walk change the key; the other five cells have no J != 0 step; no key-kept J != 0 step seen at n = 5, 6. Referee: not new (E-140, E-142), headline mostly vacuous, "refused" was a re-implementation. Response: rescoped, "cannot leave its class" deleted, library-guard cross-check 44 of 44 agree. **Accept, narrowed**. Promoted: E-293. The R-candidate "a guarded walk stays in one class" is NOT supported by this: it needs n = 7 class 1 (13 key-kept J != 0 steps, E-201) or record 6 of class 2 (E-212).

## scholar -- necessity step sourced? (referee theorist: minor revision, answered)
Claim (after rescoping): recheck of E-272 with no new result. Ladkani math/0610685 v1 Cor. 3.13 / Lemma 3.14 / Cor. 3.15 and math/0702060 v3 match as recorded; Phi-equivariance appears only inside the proof of Lemma 3.14; finite-global-dimension congruence is stated in 0702060, not ours; Nikulin Cor. 1.13.3 still UNVERIFIED (no Cartan-level NO relies on it). **Note**: no E-entry (no run); draft promoted to `research/literature/ladkani-necessity-recheck-E272.md`.

## Consequences
None for F/R. H-024 stays SUPPORTED with a stronger NO-side note (E-292). The two blocking checks of the round 108 ledger: F-candidate H-024 now has its planted/constructed NO at lam >= 1.661 but not at lam < 1.5; R-candidate still blocked (n = 7). Both go to the phase-3 conference (round 112).

## Questions for the steering committee
None.

## Decisions taken for the steering committee
Step 0.5: nothing to decide. No overnight run approved (the n = 7 class 1 walk of E-201 stays an OVERNIGHT proposal, 27 min per E-201; not added to Menu 4 this round).
