# Review of workshop/rounds/055/scholar.md

referee: skeptic · round: 055
verdict: minor revision

## Reproduction

`timeout 600 .venv/bin/python workshop/rounds/055/scholar_chz_nonmonomial.py`: 0.1 s, output identical to the submission (supp soc P_1 = {4,5}, tail-maximal ends {5}; J_1 dim 1, others 0; Prop 3.5(2) closure False, Cor 3.6(2) True). Hand check of the algebra: ab+cd is nonzero in P_1 (I only contains (ab+cd)e), killed by e, so it is a socle element at vertex 4; ab alone is not in I, so the path abe is tail-maximal ending at 5. So Cor 3.6(2) holds at S={4} while Prop 3.5(2) fails. I is admissible (length 3). Confirmed.

## True?

Checked against the LaTeX:
- AI lines 945-946 (Thm 2.31), 982-1001 (Thm 2.32), 1651-1665 (Thm 2.50): quotes are faithful. The proof printed is for (a) with (b) "dual"; the submission's reading for mu^- is the dual and is fine. Left/right is consistent: 2.32(b) is mu^- with a right D-approximation and Hom(g,D) injective; 2.50(a) says the Okuyama-Rickard complex is mu^-.
- "Silting but not tilting iff some J_i != 0": follows from 2.31 plus 2.32(b), because the proof splits Hom(T,T[!=0]) into Hom(N,D[!=0]) (exactly the injectivity), Hom(D,N[!=0]) automatic, Hom(D,D[!=0]) tilting, Hom(N,N[!=0]) automatic given the first (lines 1003-1036). So the J=0 premise is whole-T, and the "half shown" gap of E-357 is closed. Accepted.
- AI 2.50(c) vs CHZ Prop 3.5(2) at S={v}: CHZ Phi_+ = supp soc P(S) (CHZ 1269-1270), so Phi_+(S^c) in S^c iff v is not in supp soc P_j for j != v, i.e. Hom(S_v,P_j)=0. Same as 2.50(c)(ii) with e=e_v, v loopless. Accepted. This is a cleaner fact than the submission suggests and should be stated as a one-line lemma, with the loopless assumption.
- CHZ Cor 3.6 has no "monomial" (grep: no output); Ex 3.4 at line ~1298 is where tail-maximal paths = socle is asserted. Counterexample confirmed above. Caveat: "iff for every monomial kQ/I" rests on the submission's own argument (monomial implies every socle element is a path, hence tail-maximal), which I find correct but did not find printed in the paper; say so.
- Table rows spot-checked beyond the headline (EXPERIMENTS.md line numbers): row 5 (:245 text matches), row 7 (:364 "minimal left approximation of P_v is the sum of P_tb", and :419 "left add(A/P_v)-approximation"), row 8 (:364 "J_i = H^{-1} of the mutation cone"), row 10 (README.md:26 "iff for 'this mutation is a derived equivalence'"), row 12 (README.md:37 "stated for kA_n/I"), row 9 (1504 summary :133 "condition Hom(T,T[<0]) = 0 of AI Thm 2.32"), rows 13, 14 (:774, :900). All quoted text exists at those lines; corrections are right. Row 7: it is a wording/convention slip at worst, since for right modules the arrow-cover P_tb -> P_v is a right approximation (the submission already notes E-355 says "right"); not a mathematical error in the J formula.
- Edits to the two literature summaries: they are already in HEAD (d1c0005), so `git diff` on research/literature is empty; I read the r055-marked passages instead. They match the table (whole-object paragraph, N = cone[-1] in degrees 0 and 1, Thm 3.1 hereditary/canonical only, Cor 3.6 correction with the 5-vertex example). No misreading found. Minor: 2509 summary line 131 says "A tail-maximal path is always a socle element": true. "Prop 3.5 is proved for every artin algebra": matches the statement.

Not an error but a gap: row 14 says "V (mechanism)" yet the E-068 parent was not rerun; fine as stated.

## New?

`grep -rn monomial research/EXPERIMENTS.md`: E-068 (:774 area), E-124 and :900 carry the caveat from memory; E-094 (:643) concerns cord members, not CHZ. No earlier record of "the paper has no monomial hypothesis, Ex 3.4 needs it" nor of the 5-vertex example; the 2.50(c) = Prop 3.5(2) at S={v} identification is not in the record either. Genuinely new but small; mostly a rediscovery/citation audit as labelled.

## Evidenced?

Yes for the quotes (line numbers verifiable, quotes verbatim to what I read). The 5-vertex counterexample is re-runnable and hand-checkable. The only weak point is Claim 2's "automatic" steps, which are read off a proof printed only for (a); the submission says "mirrored", which is honest.

## Scope

Title says "Aihara-Iyama 2.31/2.32 and the J-criterion are verified against the LaTeX"; fine. "Cor 3.6 ... needs one" is right. The statement "for Nakayama algebras Cor 3.6 is an iff" is our derivation, not a paper claim; mark it so.

## Required for acceptance

1. Say explicitly that "Cor 3.6 is an iff for monomial I" is the submission's derivation (not printed in CHZ), one sentence.
2. Row 7: label as wording slip (right vs left approximation, module convention), not a convention error; or state which.
3. State the AI 2.50(c) = CHZ 3.5(2) at S={v} lemma with its hypotheses (v loopless, any finite dimensional algebra) as a separate numbered line so the chair can cite it.
4. Chair: the STEERING request can be cleared on the strength of this review; apply the table's research/ corrections (not done by the scholar).
