# Round 131 -- proceedings
Cycle 2, phase 3 (round 3 of 4), ordinary. Scholar, theorist, maverick; referees skeptic, experimentalist, toolsmith (all minor revision, all answered in step 3.5, checked by the chair against the review lists).

## Submissions
- **Theorist** (T10/T5, request of 130): at a gate-admitted vertex dim e_i A e_v <= 1 forces J_i = 0 (this part is E-128 L1, a rediscovery); new: with all dims <= 1 and I admissible the simple-path library gate is exact even on cyclic quivers; J != 0 at an admitted step needs dim >= 2 and a non-monomial piece of I (per row 29/29, 92/92). Referee: minor revision (false side remark on nilpotence, (d) not per row). Response did both, added E-329 bag record check. **Accept, narrowed** to what is new; promoted as **E-338**.
- **Scholar** (T3/T8, H-025 prior art): no paper states either direction; LMR Thm 6.3(3) proves N_15(3) ~ N_15(8) in Fuchsian type <2,3,10], a non-hypersurface complete intersection, and the scholar's argument (Serre functor (a)[1], Tate-Gulliksen) makes them non-FCY with Phi of order 144. Referee: minor revision (wrong claim about the record, base change, source for S, wording). Response did all. **Accept, narrowed**: the non-FCY step is an unpublished argument refereed by one persona; promoted as **E-339**, H-025 status line gets a pointer, stays SUPPORTED. Record correction: literature/2112.15587 calls Thm 6.3(3) conditional; (2,3,10) is proved outright (summary file not edited this round).
- **Maverick** (T5, R-candidate blocker): no LNA walk checked reaches a cyclic node (n = 6 d <= 5, n = 7 d <= 4, n = 8 d <= 2), so the blocker as stated cannot be met there; on random acyclic roots the library gate over-admits 39 times in 18 145 cyclic monomial rows (Hom-verified); against what E-310 (15) measured, 1 row. Referee: minor revision; response narrowed wording, fixed the 34/4/1 split, checked search.py:446, reported 6 018 swallowed `CyclicDescendantsError`. **Accept, narrowed**; promoted as **E-340**.

## Consequences
- E-310 (15)'s cyclic-parent clause ("gate decides tilting") is false for the library functions on cyclic monomial parents (E-340); R-candidate to be written at library scope at conference 132 (scope: monomial, n = 5-7, library functions, not LNA walks).
- E-337's null is explained: thin balls have no J != 0 step by E-338 (a), so the option's zero cost there is a theorem, not evidence about the option. T10 still needs a non-thin test.
- H-025 now has a proof-level (but unpublished-argument) instance at n = 15; E-245/E-246 object-level groups at n = 10-12 are not Fuchsian and are untouched.
- literature/2112.15587 summary mislabels Thm 6.3(3) as conditional: owed correction (scholar).

## Questions for the steering committee
None. Recommended defaults recorded below.

## Decisions taken for the steering committee
- Round 130 left no open question; nothing to settle (step 0.5).
- No overnight run proposed; `refuseNonTilting` default stays off.

## Next
Round 132 is the phase-3 conference: ledger names F-candidate E-323/E-332 crossings and the R-candidate (E-340); E-338 (thin => J = 0) and E-339 are the new material.
