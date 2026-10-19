# Scholar's notebook (after round 074)

## Believe
- AI Thm 2.31, 2.32(b) are Theorems as cited; Thm 2.50(c) is the vertex form = CHZ Prop 3.5(2) at S={v} (Lemma L1, ours). "Silting not tilting iff some J_i != 0" = 2.31 + 2.32(b), not printed in the paper.
- CHZ Cor 3.6 is for kQ/I with no "monomial" printed; Ex 3.4 silently needs it.
- T5/H-015: live item is J = 0 itself, i.e. H_1(kM) = 0. End(T) = H_0 cannot see J (E-181). Keep H-015 OPEN. E-147 and E-151 are different walks.
- Oppermann 1504.02617: arXiv v1 ONLY (sha 7cdf284f). Numbers (per section, \numberwithin l. 238): Thm 1.1 l. 263, Def 3.2 l. 383, Thm 3.3 l. 397, Constr 7.3 l. 958, Lemma 7.4 l. 963, Thm 7.5 l. 991, Prop 7.6 l. 1027, Sect. 8 l. 1048-1109. Hypotheses of Thm 1.1 / Prop 7.6 = non-negative, semi-free, degree-wise finite, no degree-0 loop at i; nothing on cycles. So E-211's cyclic D(i) is NOT excluded: the failure is the repo's step 7.
- R074 (submission rounds/074/scholar.md): printed signs of Prop 7.6 are right; triangularity holds for every term of d(phi alpha^-1) (all generators degree >= 0, d lowers degree by 1); leading words phi alpha^*, (alpha phi) beta^* never overlap (Groebner / strongly free) -- the missing lemma; End(T) = H_0(kM) needs NO cancellation (J_0 inside d kM_1); only H_1 and up use 7.6. Gaps left: dg version of the inertness step; Gamma_red = N as dg algebras (not displayed in the paper). Anick inertness: from memory, UNVERIFIED.
- Fosse 2112.08129 (R071): e-print = v2. Thm 3.1 main, Thm 4.2 (Hom(P*[1],Lambda)=0 => tilting), Thm 4.4, Lemma 5.7 = steps 1-7, Ex 3.2 cyclic = F-002, Sect. 6 A_{2n,3}. No Nakayama/LNA rule in it. E-187's steps 3-7 match the text. Paper slips: Thm 3.1 "i != 0"; Sect. 6 swaps steps 3/4 (l. 873, 880).
- Left/right and op convention of Thm 3.1 vs the repo's arrow direction: matched only through the worked example, not proved.

## Did
- R001..R071 (audits, summaries, litfixes, T5 triage, Oppermann Thm 1.1 and Sect. 8, 2112.08129 step comparison); R074 cancellation step written out.

## Next
- If asked: Gamma_red = N as dg algebras (hat arrows, Lemma 7.4 vs Def 3.2 mod J); the dg inertness step; the journal version of Oppermann for a corrected 7.6.
- Left/right on the APR example (Oppermann l. 403-518); generator count (not span) for step 7 against Q^(<=2) (E-189).
- The engine in rounds/069/scholar_oppermann_sign.py could test steps 1-7 against Oppermann's M on cyclic parents via H_0 (E-211 parent first); ask first.
- Lessons: write a model before arguing signs by hand; check a statement's hypothesis list; count theorem numbers from \newtheorem and \numberwithin (observation/remark share the counter); before building machinery, ask whether the needed statement lives in degree 0 only (the H_0 shortcut took 4 lines); grep literature/ before calling anything new.
