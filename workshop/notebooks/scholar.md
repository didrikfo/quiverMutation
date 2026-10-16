# Scholar's notebook (after round 069)

## Believe
- AI Thm 2.31, 2.32(b) are Theorems as cited; Thm 2.50(c) is the vertex form = CHZ Prop 3.5(2) at S={v} (Lemma L1, ours). "Silting not tilting iff some J_i != 0" = 2.31 + 2.32(b), not printed in the paper.
- CHZ Cor 3.6 is for kQ/I with no "monomial" printed; Ex 3.4 silently needs it. Iff for monomial I is our derivation.
- AI Prop 2.3 (Keller): algebraic triangulated categories. "Left approximation" for the out-arrow cover P_tb -> P_v is a right approximation (AI Def 2.30) in wording.
- T5/H-015: live item is J = 0 itself, i.e. H_1(kM) = 0. End(T) = H_0 cannot see J (E-181). Keep H-015 OPEN.
- E-147 and E-151 are different walks; no step-set equating without a key lookup.
- Oppermann 1504.02617: arXiv has v1 ONLY (e-print sha = our file, 7cdf284f); journal version UNVERIFIED. Line refs in the v1 .tex: Thm 1.1 l. 262, A = degree-0 arrows out of i l. 308, Def 3.2 l. 383-391, Thm 3.3 l. 397, Prop 7.6 l. 1027, Sect. 8 Thm l. 1048, Lemma l. 1055, Prop l. 1109.
- R069: Sect. 8 (d^2 = 0) uses only: non-negative grading, A = degree-0 non-loop arrows at i, d(alpha) = 0 by degree (l. 308). Lemma l. 1058 correct (typo t(phi) for t(p)); proof l. 1104-1105 mis-simplified ((-1)^{||p_2||} vs (-1)^{|p_2|}; they differ by -(-1)^{c_i(t p)}). Def 3.2 signs pass d^2 = 0 on 7 082 random letters (deg <= 2); flipping the last term of d(a phi b^-1) breaks it. Scripts workshop/rounds/069/scholar_oppermann_{sign,lemma}.py.
- Oppermann Thm 1.1 quasi-iso / cancellation (Prop 7.6) still a sketch (R066).
- Oppermann P_i^* = Cone(P_i -> sum P_tail(alpha)) matches AI mu^+ at the cone-triangle level; module-side convention unsettled.
- R061 (G2): step 7 spans red(Z_1 e_i) in my model; acyclic only (E-183, E-193, E-199).

## Did
- R001..R066 (audits, summaries, litfixes, T5 triage, Thm 1.1 read, version check); R069 sign/d(alpha)/d^2 test of Sect. 8 (submission rounds/069/scholar.md).

## Next
- Left/right on the APR example (l. 403-518): skeptic's; I can read for the expected quiver if asked.
- If Thm 1.1 is needed as a theorem: prove 7.6 (degree filtration) or get the journal PDF.
- The engine in scholar_oppermann_sign.py could test the repo's rewrite (steps 1-7) against Oppermann's M on cyclic parents, via M's H_0; ask first.
- Lessons: write a model before arguing signs by hand (hand reading of the printed "-" nearly mis-set my control); a statement hypothesis list plus "recall that" in Sect. 3; name theorems by line.
