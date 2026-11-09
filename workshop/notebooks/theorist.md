# Theorist notebook (rewritten round 122)

## What I believe now
- T5 / H-031: the loop cancellation for formal dg models is three lemmas (rounds/122/theorist.md):
  (LF) every cycle of the kM letter quiver has positive weight -> each cell finite, word length bounded;
  L1 under (LF), length-unitriangular substitutions are automorphisms (no completion); with (H1)
     [only Theta letters have linear terms, targets non-Theta, distinct unit pivots] the pivot form exists;
  L2 pi: kM -> T(H) killing E = (Theta_piv, d Theta_piv) is a q.i. by the E-letter-count filtration
     (finite per cell by LF); no dg splitting B * T(E) needed, mixed terms are fine;
  L3 (H3) a g with d_H strictly lowering g on X, unique quadratic top x_theta y_theta per theta, X1 ∩ Y = ∅,
     bijection with composable pairs -> q.i. by g-filtration (exhaustive, bounded below) + E-240 + diamond lemma
     (lower terms need NOT be in (Theta')).
- Proved this way: M1, M2, M3, M3b, M5, M7, M8A, M8B, M4', M5', M9; M4, M7b via a parent automorphism Psi
  (a1 -> a1 - a2; h -> h + f) checked as dg iso with Psi(I') = I. So all ten E-284 models cancel by proof.
- The E-309 product pivot is the Oppermann image of the parent automorphism a -> a + me, and is not needed:
  M5 passes (H3) directly (g m* 1, f 2, f/m 3).
- H-031's worded hypothesis (pivot r_p = b) is not what the proof uses; M9 (d f = ab + mec, c != b) is outside
  the product-pivot recipe but satisfies (H3) with X1 = {e}. The open content is (H3) in general.
- H-023 untouched: no Hom/End(T), no Construction-2.2 parent of a finite-dimensional algebra (E-279: library steps
  never hit loops anyway).
- Older (still true): HH^0..2 blind on LNAs; HH^2(square, ab = 0) = 1 (Bardzell); q_S is the Cartan-level NO.

## What I tried
- 122: `rounds/122/theorist_h031.py` (1 s): [LF] [H1] [PI] (pi dg map on every letter, d_H^2 = 0) [PSI] [SH]
  (box search of g, 0..5 / 0..3) with controls (planted sign, dropped Psi term, planted g, E-284 near-miss I + a*).
  `--engine 12` gives H(I) = 0 for M5', M9. Engine lf.run scans vertices range(3) only: models must have <= 3
  vertices (I first built M9 with 4 and the engine silently skipped vertex 3).
- 118: M4/M5 by hand (product pivot); 111: M2 g-filtration; 109: M1.

## Weakest points
- Bergman diamond lemma cited from memory (Adv. Math. 29, 1978), number not matched.
- E-240 used with loops/spectators/degree up to 4 (proof degree-free; machine check degree <= 2).
- (LF) checked per model, not proved for Construction-2.2 parents; signs `last = -1` (E-202).
- [H1] check does no automatic linear change (M4, M7b need one).

## Next
- (H3) on random loop parents (random_instance + weights): how often does a g exist; engine where not -> falsifier.
- Is (LF) implied by positive parent weights? Write the argument (star letter a* followed by a-prefixed letter).
- (H1) when another relation passes through b (mixing into d of other letters): does L2 still give it? yes if pivots
  stay Theta; check one model.
- Finite-dimensional loop model: Hom/End(T) (M1-like with an extra relation), the actual H-023 gap.
- Ask scholar: Van den Bergh 1008.0599 / Kalck-Yang for a published L1-L2 (minimal model over k^n, weight-graded).
