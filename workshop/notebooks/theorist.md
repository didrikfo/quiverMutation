# Theorist notebook (rewritten round 131)

## What I believe now
- Thinness (round 131): gate + d_i = dim e_i A e_v <= 1 => J_i = 0. This is E-128 L1 (my own, round 034): gate => dim J_i <= d_i - 1.
  - Hypotheses: paths span e_iAe_v (I inside R), the gate quantifies over ALL nonzero paths into v, out-arrows counted
    as arrows, v has an out-arrow. Not used: acyclic, monomial, field, out-degree, key guard.
  - New: if all d_w <= 1 (w != v) and R is nilpotent mod I, every nonzero path into v is vertex-simple (first repeated
    vertex w gives p2, C p2 independent by nilpotence), so the library's simple-path gate is exact at thin vertices, cyclic too.
  - Exception shape: gate-admitted J != 0 needs d_i >= 2 AND e_i I e_t(beta) not spanned by paths; J generator has support >= 2
    (E-080 / E-318 shape p - lambda q). Monomial: gate <=> J = 0 at any d (E-208, E-310).
  - E-337's "isTiltingAtVertex never False" is implied by "d <= 1 everywhere", not a second observation. 13 F-041 balls,
    6 070 steps: all on acyclic quivers.
- T10 mechanism: "J != 0 <=> non-tilting irreducible silting mutation" is a theorem on record (E-130, E-171, E-138; my check E-333).
  It does NOT give child = End(T) (that is H-023), and non-tilting does not imply out of class (E-316's 19).
- E-328's any-slot certificate needs a finite list of indecomposables (type D only): useless at n = 11, 12.
- T5 / H-031 route: (LF) + (H1w) + (H3) => kM -> N q.i. (LF) proved under (C); (H3) open in general.
  Monomial cyclic Nakayama (E-326) all hold; g is not path length; (C)/(LF) depend on the presentation.
- Older: (C) whenever |A| = 1 or the A-weights are equal; (LF) <=> potential exists; HH^0..2 blind on LNAs.

## What I tried
- 131: `rounds/131/theorist_thin.py`: an independent exact engine for acyclic kQ/I (ideal part from u r w, row-reduced over Q),
  d_i, gate, J_i; sweep against the library (0 disagreements on 3 200 vertices), planted failures, cyclic monomial simplicity check.
  `theorist_thincyc.py`: acyclicity counter on F-041 balls via round-130 toolsmith script (~5 s per pair).
- 129: `rounds/129/theorist_jsilt.py`: two-term complexes, Hom(T,T[m]), End^0(T) Cartan; reusable `homdim`.
- Lesson (third time): grep research for the claim BEFORE proving. "d_i - 1", "J_i", "silting", "Hom(T,T[-1])".
  This round the assigned question was E-128's title.
- 126: `theorist_cycparent.py`, `theorist_creach.py`. 123: `theorist_lnaparent.py`, `theorist_shift.py`.

## Weakest points
- Child = End(T): Cartan level only (mine, E-161), quiver level on 13 edges (E-167). This is the real gap behind H-023 and T10.
- Library `isInIdeal` on cyclic quivers (E-329): my (c) is mathematical; the library could still misjudge membership.
- Minimality of the approximation with parallel arrows: my own argument. (H3) analytic half; Bergman from memory.

## Next
- H-023 at algebra level: the End(T) multiplication table from the Hom^0 basis (engine of 129), then an algebra isomorphism test
  against the child on small n, for J = 0 steps first and then J != 0.
- Where do d >= 2 steps occur on walks? (E-332 at depth 6, E-316 mixed). Is there a structural reason F-only depth-3 balls stay thin?
- T7 / H-010 end-only gap-2 argument (still owed). T5: (H3) on a non-monomial or |A| >= 2 parent.
