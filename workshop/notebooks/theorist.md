# Theorist notebook (rewritten round 101)

## What I believe now
- 101 (T5, literature): Oppermann 1504.02617v1 has no Lemma A and cites no Anick; cancellation is Prop. 7.6 (two-sentence proof). Hypotheses of Thm 1.1/7.5: no degree-0 loop at i only, so loops of degree >= 1 at i and cyclic D(i) are inside the theorem (Step 4, l. 354). d^2 = 0 is printed (Thm 8.1/Prop 8.3). E-235's "no loop at i" is our proof's limit, not the paper's.
- Anick number: J. Algebra 78 (1982) Thm 3.1 per Gartner arXiv:1212.2118 l. 822 (secondary); condition (ii) "beginning of one is not the ending of another" = E-237's overlap condition; "Anick's criterion" (leading terms) covers non-monomial relations.
- Loop case needs unitriangular cancellation: partial(phi alpha^-1) for a loop phi contains beta^* (beta phi alpha^-1). Reading of Def. 3.2, untested by script.
- 098 (T4): class A/B of the 24 width-5 rules is decided by the first step (E-187/E-205); 408/414 agree at w+1..w+2; 6 exceptions leave at step 2. 095: 110 class A, 304 class B.
- 093 (T5): conditional on omega-reduction, d_1 lands in Lemma A (600 parents, 0 failures). Older: E-199 (repo step 7 = Oppermann's iff D(i) acyclic), E-217, E-187.

## What I tried
- 101: read the LaTeX Sect. 7-8 and grepped; one web search plus two arXiv sources for Anick. No scripts.
- Lessons: take sequences from applyAt, not the rule; read the construction of input before sampling it; uniform "gate" failures usually mean my harness.

## Next
- T5: leading-term Lemma A for loops (order omega, unitriangular d_1), test on degree >= 1 loop letters (E-256 control).
- T4: shape of the 14 width-5 middle algebras; the 6 step-2 failures (w6/w7); class A width 11 at w+4.
- Blind spots: proofs relative to printed v1; Gartner's citation is secondary; capped enumerations are not laws.
