# Scholar's notebook (after round 063)

## Believe
- AI Thm 2.31, 2.32(b) are Theorems as cited; Thm 2.50(c) is the vertex form = CHZ Prop 3.5(2) at S={v} (Lemma L1, ours). "Silting not tilting iff some J_i != 0" = 2.31 + 2.32(b), not printed in the paper.
- CHZ Cor 3.6 is for kQ/I with no "monomial" printed; Ex 3.4 silently needs it. Iff for monomial I is our derivation.
- AI Prop 2.3 (Keller): algebraic triangulated categories. "Left approximation" for the out-arrow cover P_tb -> P_v is a right approximation (AI Def 2.30) in wording.
- T5/H-015: live item is J = 0 itself, i.e. H_1(kM) = 0. End(T) = H_0 cannot see J (E-181). Keep H-015 OPEN.
- E-147 and E-151 are different walks; no step-set equating without a key lookup.
- Oppermann 1504.02617: one file mutation.tex, 1173 lines, sha256 7cdf284f..; \numberwithin{theorem}{section}, so main theorem = 1.1 (Sect. 1). Sect. 7 items counted by hand, name them by line (895-1043).
- R063 (referee-corrected): Thm 1.1 hypotheses invoked in Sect. 7: only no degree-0 loop at i (lines 264, 321); Step 4 (l.354) = loops of degree>=1, longer cycles via Steps 2/3; step 5 (1027-1042) is a sketch, Sect. 8 unread. Hypothesis list duplicates lit file:36-40 and E-181; new = proof read + Q^(<=2) for H_0 (by hand). Left vs right unsettled.
- R061 (G2): step 7 spans red(Z_1 e_i) in my model; 35 + 573 random pairs equal (rounds/061/scholar_g2.py); acyclic, span level only.
- Ex. 4.2 mutates at vertex 2 where step 7 is empty; tests steps 4-6 only. Ex at line 403 (APR tilt at 1) is the one to run for the convention.

## Did
- R001..R058 (citation audit, summaries, counterexample, litfixes, T5 triage); R061 G2 script; R063 read the proof of Thm 1.1 (Sect. 7), submission rounds/063/scholar.md. No script.

## Next
- Ask whether the chair put mutation.tex in sources/ (version number unread; get it from the arXiv listing).
- Read Sect. 8 (well-definedness, line 1044+) only if someone doubts d~^2 = 0.
- Left/right test on the APR example with the repo (hand to skeptic); generator count and parallel arrows for G2 (theorist lemma).
- Lessons: the hypothesis list of a theorem is in its statement plus one "recall that" in Sect. 3; check where each is used in the proof. Name theorems by line when numbers are only counted.
