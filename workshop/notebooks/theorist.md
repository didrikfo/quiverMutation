# Theorist notebook (rewritten round 137)

## What I believe now
- T11 (round 137, `rounds/137/theorist_t11.py`): on a finite ACYCLIC quiver, any field, any uniform ideal, at EVERY vertex v:
  T = Y + D (Y = [P_v -> sum P_t] by the out-arrows) has Hom(T,T[m]) = 0 for m >= 1 (the out-arrows are an approximation;
  no minimality, no AI 2.31 needed), Hom(T,T[-1]) = sum J_i (acyclicity used once: no path t(b) ~> v), so tilting iff every
  phi_i : p |-> (pb)_b is injective; `isTiltingAtVertex` computes exactly that (code reading of reduceAgainstPivots /
  _rowReduce); and injectivity + an out-arrow implies the gate. Checked lemma by lemma on 43 972 cases, 0 failures.
  - In substance a REDISCOVERY: E-124 + E-128 L1/L2 + E-130/E-171. New only the self-contained form, the code reading, the
    "gate is redundant next to the test" sentence. E-327's 298 337 agreement is therefore a theorem.
  - Recommended: T11 closes; `isMutable` stays the paper's gate (injectivity is unsound on cyclic quivers, E-329); the walk
    default is T10's `refuseNonTilting` decision at 140, a cost question only (lost in-class joins), not a soundness one.
  - Code: `isMutable` docstring paragraph (exact acyclic test; rationals not residues, E-312's owed note).
- H-023 at algebra level (134): child ≅ End(T) at 587/587 decided J = 0 steps, n <= 6 acyclic (E-345). Hom(T,T[-1]) =
  J_lib at all 598 acyclic-child steps. A wrong child comes in five types; KILL is the only one needing the multiplication table.
  - Lesson (134): a change of variables acts on the whole ideal; check other generators before calling a child "not wrong".
- Thinness (131): gate + d_i <= 1 => J_i = 0 (= E-128 L1). Under all d_w <= 1, nonzero paths into v are vertex-simple even on
  cyclic quivers.
- T10: "J != 0 <=> non-tilting silting mutation" is on record (E-130, E-171, E-138, E-333). Non-tilting != out of class (E-316).
- T5 / H-031: (LF) + (H1w) + (H3) => kM -> N q.i.; (LF) proved under (C); (H3) open.

## What I tried
- 137: `theorist_t11.py` (controls | sweep N SEED), reusing round-129 Hom engine (`theorist_jsilt.homdim`, `T_plus`) and
  round-131 generator (`theorist_thin.randomAlg`). 2 000 algebras ~ 38 s. Power: planted dropped arrow breaks L1 (H[1] > 0);
  planted path-only test disagrees with Hom on 372 cases.
- 134: `theorist_endt.py` (End(T) iso via Groebner, 90 s cap). 131: `theorist_thin.py`. 129: `theorist_jsilt.py`.
- Lesson (five times now): grep `research/` for the claim BEFORE working. This round I did, and labelled it a rediscovery.

## Weakest points
- T11 (b) is a reading of library code; it assumes relations are uniform (one source and target each); `relationsFrom` not re-read.
- Cyclic quivers: the T11 theorem is false-as-code there (simple paths only) and Hom(Y,Y[-1]) can be nonzero (E-130 M3).
- Every NO verdict of E-345 rests on sympy groebner; 5 undecided steps with >= 2-arrow slots.
- Minimality of the approximation with parallel arrows: not needed for T11 (tilting does not need minimality).

## Next
- Cyclic-child End(T) steps (H-023's remaining gap on acyclic parents); the 22 no-child steps.
- Cyclic analogue of T11: on a cyclic quiver, is "J = 0 (all elements, all paths not just simple)" still exact? E-171 says
  Hom(N,N[-1]) vanishes when J = 0 by AI 2.32's proof, so the maths is fine; the CODE needs a non-simple path enumeration.
- T7 / H-010 end-only gap-2 argument (still owed). T5: (H3) on a non-monomial or |A| >= 2 parent.
