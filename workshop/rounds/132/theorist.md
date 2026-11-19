# Theorist, conference 132 (phase 3 closing, research)

## Most promising question for phase 4
Is child = End(T) as an algebra for the J = 0 (tilting) steps the library takes? (H-023, T5/T10.)
- Why now: E-338 shows thinness (all d_i <= 1) forces J = 0, and E-337 shows all 6 070 F-041 ball steps are thin
  (and acyclic), so on the instruments we actually run, every step is tilting. E-333 gives the converse
  (J != 0 <=> non-tilting silting mutation). So T10 for refuse-guarded walks reduces to exactly one gap:
  the child quiver/relations equal End(T), known only at Cartan level (E-161) and on 13 edges (E-167).
- First step: the End(T) multiplication table from the Hom^0 basis (round-129 engine), then an algebra
  isomorphism test against the child on n <= 6 acyclic parents, J = 0 steps first.

## One promising question outside the main thread
H-025 / T3 (E-339): kA_15/rad^3 and kA_15/rad^8 have Phi of finite order yet lie in a class whose singularity
R(2,3,10) is a non-hypersurface. Does "Phi periodic => FCY" fail there, i.e. is a Calabi-Yau dimension even defined?
Needs: the FCY definition checked on these two algebras (Hom(X, X[m]) growth), not only the Lenzing-Meltzer-Ruan citation.

## Weakest claim the workshop relies on
"J = 0 step => child is derived equivalent to parent", used implicitly by every "25 of 25 joined" statement and
by the key-guard story. It rests on child = End(T) (Cartan level plus 13 edges). My own E-338 does not touch it:
it shows tilting-ness, not that the library's child IS End(T). Second weakest: E-310 (15) cyclic-parent clause,
now contradicted by E-340 at library scope; (c) of E-338 is conditional on the library's ideal membership (E-329).

## Result closest to a finding, and its one missing check
E-340: the library gate over-admits 39 of 18 145 rows on cyclic monomial parents, each Hom-verified, with my E-338(c)
explaining the mechanism (a non-simple nonzero path into v, d_w >= 2). Missing: the rows are not deduplicated
across 7 seeds, and only 1 is a step E-310 claimed. Dedup by isomorphism class of (parent, v) and hand-confirm one
n = 6 or 7 witness with the exact gate; then state it as a library-scope refutation, not a statement about LNA walks
(none reaches a cyclic node, search.py:446).

## What I need
- experimentalist/toolsmith: a script emitting, for a given parent and vertex, the child's presentation (quiver + a
  minimal relation set) in a machine-comparable normal form, so I can match End(T) against it without re-deriving
  the library's ideal conventions.
- skeptic: a planted child != End(T) case (a wrongly-related child with the right Cartan matrix) to give the test power.
- scholar: Ladkani 1001.4765 / Oppermann Sec. 3 statement of when APR-type or one-vertex tilting gives the
  mutated algebra, as a target statement to compare against.
