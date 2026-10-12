# Scholar's notebook (after round 061)

## Believe
- AI Thm 2.31, 2.32(b) are Theorems as cited; Thm 2.50(c) is the vertex form = CHZ Prop 3.5(2) at S={v} (Lemma L1, ours). "Silting not tilting iff some J_i != 0" = 2.31 + 2.32(b), not printed in the paper.
- CHZ Cor 3.6 is for kQ/I with no "monomial" printed; Ex 3.4 silently needs it; 5-vertex counterexample (round 055 script). Iff for monomial I is our derivation.
- AI Prop 2.3 (Keller) is for algebraic triangulated categories; Rickard summary already says so.
- Wording: "left approximation" for the out-arrow cover P_tb -> P_v is a right approximation (AI Def 2.30).
- T5/H-015: the live item is the premise J = 0, not the identification; End(T) = H_0 cannot see J (E-181). Keep H-015 OPEN.
- E-147 and E-151 are different walks; do not equate step sets without a key lookup.
- R061 (G2): repo step 7 spans red(Z_1 e_i) of my model of Oppermann's degree-2 arrows (span level; sign/transpose only read from the paper, review R061). Steps 4,5,6 read off Def 3.2 and the Sect. 7 proposition (phi alpha* = -(-1)^|phi| (d phi)/alpha), same signs, not computed. Computed: 35 hand pairs + 573 random pairs, ALL EQUAL (`rounds/061/scholar_g2.py`).
- Oppermann Ex. 4.2 mutates at vertex 2 where step 7 is empty; it tests steps 4-6. The degree-2 arrow x is at vertex 1.
- Oppermann LaTeX = one file mutation.tex, 78 817 bytes, e-print 1504.02617. Theorem numbers counted from \newtheorem (shared counter by section): main theorem 1.1, Def 3.2 (differential), Ex 4.1 and 4.2, proposition in Sect. 7. Not compiled.

## Did
- R001..R058 as before (citation audit, summaries, counterexample script, litfixes, T5 triage).
- R061: fetched LaTeX, read Sect. 3, 4, 7 (proof), derived the step-7 identification by hand, wrote the script, submission `rounds/061/scholar.md`.

## Next
- Generator count vs minimal count (review item 5). Ask theorist for the one lemma I did not prove (minimal degree-2 generators = Z_1 e_i mod Z_1 A and left multiples).
- Test G2 on real walk steps (n=7 LNAs, E-151 failing steps) and parallel arrows; needs arrow-keyed words.
- Check the chair applied the R058 patches; line numbers drift, copy old text.
- Read Oppermann's proof of Thm 1.1 (Sect. 7) if the T5 conditional is to be closed; so far only statement matched.
- Lessons: Ex. 4.2 looked like the test but is empty at its vertex; check which vertex an example mutates before calling it a test. EXPERIMENTS entries renumbered at merge; verify E-numbers against the file.
