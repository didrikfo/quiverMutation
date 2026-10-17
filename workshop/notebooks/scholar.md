# Scholar's notebook (after round 071)

## Believe
- AI Thm 2.31, 2.32(b) are Theorems as cited; Thm 2.50(c) is the vertex form = CHZ Prop 3.5(2) at S={v} (Lemma L1, ours). "Silting not tilting iff some J_i != 0" = 2.31 + 2.32(b), not printed in the paper.
- CHZ Cor 3.6 is for kQ/I with no "monomial" printed; Ex 3.4 silently needs it.
- T5/H-015: live item is J = 0 itself, i.e. H_1(kM) = 0. End(T) = H_0 cannot see J (E-181). Keep H-015 OPEN. E-147 and E-151 are different walks.
- Oppermann 1504.02617: arXiv v1 ONLY (sha 7cdf284f). Line refs in the v1 .tex: Thm 1.1 l. 262, Def 3.2 l. 383-391, Thm 3.3 l. 397, Prop 7.6 l. 1027, Sect. 8 l. 1048-1109. Sect. 8 (d^2 = 0) uses only non-negative grading + A = degree-0 non-loop arrows at i; Lemma l. 1058 correct, proof l. 1105 mis-simplified (E-202). Thm 1.1 quasi-iso / cancellation (Prop 7.6) still a sketch.
- Fosse 2112.08129 (R071): e-print = v2 (arXiv lists v1, v2; sha256 of .src d5b2f668...). Numbers (counted by hand): Thm 3.1 main, Thm 4.2 hypothesis Hom(P*[1],Lambda)=0 => tilting, Thm 4.4 the two impossible cases, Lemma 5.7 = steps 1-7, Ex 3.1/3.2 (3.2 cyclic = F-002), Sect. 6 A_{2n,3}. NO Nakayama/LNA rule in it; the one-step rule is ours. E-187's steps 3-7 match the printed text (only a composition-order wobble in prose). Hom(P*[1],Lambda) = {g : g.[alpha]=0} is derived by me, equals Thm 4.4 case 2 with elements for paths; for monomial parents gate is exact. Paper slips: Thm 3.1 "i != 0"; Sect. 6 swaps steps 3/4 (l. 873, 880).
- Left/right and op convention of Thm 3.1 vs the repo's arrow direction: matched only through the worked example (scholar_paper_example.py), not proved.

## Did
- R001..R069 (audits, summaries, litfixes, T5 triage, Oppermann Thm 1.1 and Sect. 8); R071 2112.08129 step comparison (submission rounds/071/scholar.md), draft summary patch for the chair.

## Next
- If asked: left/right on the APR example (Oppermann l. 403-518); generator count (not span) for step 7 against Q^(<=2) (E-189).
- If Oppermann Thm 1.1 is needed as a theorem: prove 7.6 or get the journal PDF.
- The engine in rounds/069/scholar_oppermann_sign.py could test steps 1-7 against Oppermann's M on cyclic parents via H_0; ask first.
- Lessons: write a model before arguing signs by hand; check a statement's hypothesis list; count theorem numbers from \newtheorem (observation/remark share the theorem counter, they shifted my first count); a derivation from a cone complex took 5 lines and showed the "gate vs tilting" gap is the known 1001.4765 one: grep literature/ before calling it new.
