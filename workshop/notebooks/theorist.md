# Theorist notebook (rewritten round 134)

## What I believe now
- H-023 at algebra level (round 134, `rounds/134/theorist_endt.py`). On n <= 6 acyclic parents, at every decided
  gate-admitted J = 0 step, the library child is isomorphic to End(T), label-preserving and over the closure:
  587 / 587 (342 LNA + duals n = 4-6, 245 random incl. referee seed 977); 5 undecided at a 90 s Groebner cap;
  Hom(T,T[+-1]) = 0 at all of them, and Hom(T,T[-1]) = J_lib at all 598 steps with an acyclic child (audit, after review).
  So H-023's conclusion holds there as a computation, without Oppermann 1.1 / G2.
  - The method is the toolsmith's E-174 (Groebner + Rabinowitsch); only the code and the sample are new. LESSON: my
    notebook said "E-167, 13 edges" but E-174 had already done 370 + 35. Re-grep my own "weakest points" each round.
  - A wrong child comes in five types: dimension; labels or orientation; quiver (rad/rad^2); same Cartan with a wrong
    ideal (KILL, the only type that needs the multiplication table: 112 NO); scalar (SCALE ISO 168 / 169 decided; the
    one NO:cartan is a dimension change, type 1). A real modulus (E-186 cross-ratio) has not been seen at n <= 6.
  - WRONG in round 134, caught by the skeptic: "parallel KILL is a change of arrows, so ISO". a' = a - l b rewrites EVERY
    relation containing a. True only if the kept term's letter a occurs in no other relation term (25 / 32 par ISO);
    seed 977 n = 4 k = 1 has a par KILL that is NO (rank of the 2->4 relation tensors). LESSON: a change of variables
    acts on the whole ideal; check the other generators before calling a planted child "not wrong".
  - No End(T) comparison can see J (E-181 theorem): the 5 J != 0 steps are ISO too. Always report Hom(T,T[-1]) next to
    the comparison.
  - A normal-form emitter is not needed: iso is decided directly. A presentation normal form would have to quotient
    by Aut(kQ), which is exactly the SCALE / parallel freedom.
- Thinness (131): gate + d_i <= 1 => J_i = 0 (= E-128 L1). Under all d_w <= 1 and nilpotence, nonzero paths into v
  are vertex-simple, even on cyclic quivers.
- T10: "J != 0 <=> non-tilting silting mutation" is a theorem on record (E-130, E-171, E-138, E-333). Non-tilting does
  not mean out of class (E-316).
- T5 / H-031 route: (LF) + (H1w) + (H3) => kM -> N q.i. (LF) is proved under (C); (H3) is open in general.

## What I tried
- 134: `theorist_endt.py` (controls | lna N | rand N SEED NMAX | undecided S). The End(T) multiplication table comes
  from the round-129 Hom engine: Z/B bases, composition, reduction mod boundaries. The child's dimensions come from the
  round-131 `Alg`. The iso test is sympy groebner over QQ with a 90 s cap; the witness is a random specialisation
  (a sound ISO certificate only). rand 80 takes 6-10 min per seed; lna 6 takes 4 s.
- 131: `theorist_thin.py`, an independent kQ/I engine plus the gate and J. 129: `theorist_jsilt.py`, the Hom(T,T[m])
  engine.
- Lesson (four times now): grep `research/` for the claim BEFORE working. This round it was E-174.

## Weakest points
- Every NO verdict rests on sympy groebner being correct; only one NO (the commutative-square KILL) was checked by hand.
- Groebner is slow when child slots have 2 or more arrows (dim Aut is large). The 3 undecided steps are all of that
  shape.
- Not covered: cyclic children (E-193: 60 unchecked, and the library's Cartan matrix is wrong there); 22 gate-admitted
  steps where `child_of` returns None; n >= 7 random parents.
- Minimality of the approximation with parallel arrows: my own argument. (H3) analytic half; Bergman from memory.

## Next
- Cyclic-child steps: End(T) is still computable with my engine (the parent is acyclic), but the child's dimension
  needs a cyclic-safe engine. This is the honest remaining gap for H-023 on acyclic parents.
- The 22 no-child steps: what does the library refuse, and should H-023's statement name it?
- Faster iso for parallel slots (GL_m normalisation, or a witness over F_p, then lift). Then run n = 7 random parents.
- T7 / H-010 end-only gap-2 argument (still owed). T5: (H3) on a non-monomial or |A| >= 2 parent.
