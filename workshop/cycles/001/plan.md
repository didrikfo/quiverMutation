# Cycle 001 plan (rounds 101-120)

Chair's merge of the six persona plans (`plan/<persona>.md`). Not binding; each
phase conference revises it. Written 2026-10-30. The human can steer by editing
this file or `STEERING.md`.

## Questions of the cycle
1. **Is the Cartan-level NO valid? (T3/T8 · H-024, new)** Is "derived equivalent => Phi-equivariantly Z-congruent" a published theorem, and does the E-269 test pass on certified-equivalent pairs? Five of six plans lead with this. One failure voids the 9 NOs of E-265/E-269.
2. **Is every gate-admitted J = 0 step a derived equivalence? (T5/T10 · H-023, new; H-015)** Close E-235 step 4 against Lemma A's hypotheses and the loops-of-degree->=1 case, or get a per-edge certificate that does not use the move set (70 failing n = 7 edges). No more join searches.
3. **Controls for H-010 (T7)**: a stepped-chain control (8, 10, 12) that can lower overlap, then a decision.
4. **H-017 (T6)**: defect per key class on saved n = 9 data.
5. **H-020 (T4)**: one-step rule for the 14 width-5 B rules, the 6 step-2 failures held out.
6. **H-025 (Phi periodic, not FCY)** and **H-026 (S-1 K-threshold)**: literature check only unless a slot is free.

## Phase 1 reading list (rounds 101-104)
- **Necessity of Z-congruence** (checks E-265, E-267, H-024; Nikulin 1.13.3/1.14.2 UNVERIFIED): Happel, Ladkani, Barot-de la Pena, Lenzing-de la Pena, Bayer-Fluckiger, SPLAG ch. 15. Scholar, maverick. Either a matched theorem number or "no source, proof is ours".
- **Oppermann 1504.02617** steps 1-7 and cyclic D(i), from the LaTeX (checks E-199, E-211, H-023); derived equivalence with loops/cycles (Rickard, Keller). Theorist, scholar.
- **Anick Lemma A number** (E-237, UNVERIFIED), Chan-Koenig-Liu, Aihara-Iyama on non-tilting steps (E-235 step 4, H-015). Scholar, toolsmith.
- **Independent inequivalence invariants** for Nakayama/gentle algebras (Hochschild cohomology, Avella-Alaminos-Geiss, Bessenrodt-Holm; 2312.14699 on file): a second NO test for E-259's eight pairs. Skeptic, scholar.
- **Explicit equivalent pairs at n = 9..12** (Happel-Seidel, Bobinski-Geiss-Skowronski, Ladkani): independent positives for the power control. Experimentalist, scholar.
- **Phi periodic => FCY** converse (1310.1557, 2012.11927, 2512.19873 on file; H-025), orbit counts for B_n (E-247, H-026), Nakayama vertex deletion (S-1). Scholar, theorist.
- Toolsmith: reproduction examples and any per-step tilting test in implementations (QPA, Keller applet).
- One phase-1 slot (not the scholar's) goes to the signature/genus of S for the 140 group from saved data (maverick/scholar).

## Loose agenda, phases 2-4
1. **Power control for H-024** (experimentalist + toolsmith harness, maverick, skeptic pre-registers): E-269 test on the 1 002 n = 9 same-name pairs of E-238 and 6 n = 12 within-quipu pairs. Then exact integer bound for E-269 (maverick).
2. **H-023** (theorist, toolsmith, skeptic): proof of the two joints, `checkJZero` per-edge certificate on the 70 edges, silting-engine check of the key-kept J != 0 edges (skeptic).
3. **Cheap second invariant** (maverick): Z[Phi]-module class vs q_S; a split of an E-117/F-010 class is a bug.
4. Breadth: H-017 census (experimentalist), T4 width-5 rule (theorist), T7 stepped chain (experimentalist). S-1 only with a free slot.

## Expected backlog (phase 5)
`mutationSearchDepthFirst` docstring (l. 322), `canonicalKey` and guard-inert-on-cyclic docstrings; five direct-caller script reruns; E-269 float bound if not made exact; Nikulin 1.13.3 and Anick Lemma A numbers; E-235 step 4; H-020 H3/H4; S-1 item 3; undecided overnight proposals (E-204, E-206/E-210, E-253, E-260).

## What would count as progress
- H-024 decided: sourced necessity step plus 1 002 + 6 controls, or one counterexample that voids the NOs (equal progress). F-entry or R-entry.
- H-023 reduced to stated lemmas with proofs checked on the graded parents, or one counterexample; plus a per-edge certificate/refusal for each of the 70 edges.
- Every UNVERIFIED citation in STATE either matched to its theorem number or struck.
- Library docstrings fixed and scripts rerun.

## Where the personas disagree
- **Order of work:** toolsmith and skeptic want the power control before any new Cartan invariant; maverick and scholar also propose new invariants (Z[Phi]-class, Hochschild/AAG) in the same cycle.
- **H-017 census:** skeptic would put it after the power control; experimentalist and maverick want it as cheap breadth.
- **T10:** all six reject more join searches; theorist wants a proof route, toolsmith a library certificate, skeptic an independent silting check.
- Scholar and experimentalist plans exceed the 40-line limit.
