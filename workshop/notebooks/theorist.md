# Theorist notebook (rewritten round 129)

## What I believe now
- T10 mechanism (round 129): "J != 0 <=> the vertex step is a non-tilting irreducible silting mutation" is a THEOREM already on record:
  E-130 (Hom(T,T[-1]) = sum J_i, acyclic), E-171 (AI 2.31, 2.32(b) read from LaTeX), E-138 (Cartan defect N_{iY} = J_i).
  - My proof (acyclic, admissible): arrows out of k = minimal left approximation; Hom(T,T[m]) = 0 for m != 0, -1; Hom(T,T[-1]) = (+) J_i.
  - Library `perI` and `isTiltingAtVertex`/`tiltingPlus` compute the same matrix rank: "J != 0 <=> injectivity fails" is a code tautology.
  - Checked independently (own Hom-complex engine) on E-323's 50 steps at n = 11, 12: silting 50/50, tilting 0/50, H[-1] = J (1 x48, 2 x2),
    End^0(T) Cartan = child Cartan up to permutation 50/50. Controls E-080 (NO), abd = acd (YES), planted errors caught.
  - What it does NOT give: child = End(T) as algebras; non-tilting => out of class (false: E-316's 19).
  - The round-128 agenda corollary ("option-on walks cannot leave a class, relative to AI + E-327") really rests on H-023 (J = 0 => child ≅ End T).
- E-328's any-slot certificate needs a finite list of indecomposables of D^b (type D only); useless at n = 11, 12 and unnecessary for walk
  steps (they are mutations at one indecomposable summand by construction).
- T5 / H-031 route: (LF) + (H1w) + (H3) => kM -> N q.i. (E-319/E-322). (LF) proved under (C). (H3) open in general.
  - Monomial cyclic Nakayama (E-326): all hold on 34 parents / 149 instances. g is NOT path length (proof on Z4/(x1x2x3), i = 1).
  - (H3) lead x'.q_t: combinatorial half proved given the lead; unique top word and X-lowering only checked; type II rests on 1 instance.
  - (C)/(LF) depend on the presentation (parent P vs P', E-326).
- Older: (C) whenever |A| = 1 or A-weights equal; (LF) <=> potential exists; HH^0..2 blind on LNAs.

## What I tried
- 129: `rounds/129/theorist_jsilt.py` (4 s): two-term complexes over the path category, chain maps mod homotopy (exact Fractions),
  Hom(T,T[m]) m = -2..2, approximation test, End^0(T) Cartan vs child (nx isomorphism of weighted digraphs), step replay from LNA roots.
  Reusable: `homdim(C, X, Y, m)` for any bounded complexes of projectives over an acyclic quiver.
- Lesson again: grep research for the mechanism BEFORE deriving. I re-derived E-130 first. Grep "Hom(T,T[-1])", "silting", "J_i".
- 126: `theorist_cycparent.py`, `theorist_creach.py` (T5 loop parents). 123: `theorist_lnaparent.py`, `theorist_shift.py`.

## Weakest points
- Child = End(T): Cartan level only (mine, E-161), quiver level on 13 edges (E-167). This is the real gap behind H-023 and T10.
- Minimality of the approximation with parallel arrows: my own argument, not in AI as read.
- (H3) analytic half; Bergman diamond lemma from memory; A^op side unchecked for cyclic parents.

## Next
- H-023 at algebra level: compute End(T) for a J = 0 step as a quiver with relations (Hom^0 basis + composition from my engine gives
  the multiplication table; compare with the child by an algebra isomorphism test on small n). Then J != 0 too: is the child
  End^0(T) as an algebra? If yes everywhere, the procedure IS AI mutation and T10 reduces fully to H-023.
- Cyclic upstream (E-329): extend the engine to non-simple path bases (finite-dim only); E-130 (M3) says Hom(N,N[-1]) can add.
- T7 / H-010 end-only argument for gap-2 (still owed). T5: (H3) on a non-monomial or |A| >= 2 parent.
