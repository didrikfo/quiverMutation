# Theorist notebook (rewritten round 123)

## What I believe now
- T5 / H-031 route (E-319, corrected in round 123): (LF) + (H1w) + (H3) => kM -> N q.i.
  (LF) every cycle of the kM letter quiver has positive weight;
  (H1w) every Theta letter with a linear term has a pivot (c != 0, non-Theta, distinct, in no other THETA's linear part).
     Non-Theta letters MAY have linear terms: Lemmas 1-2 never used E-319's clause "only Theta letters are linear";
  (H3) g with d_H strictly lowering g on X (linear terms allowed), one quadratic top x y per Theta', X1 ∩ Y = ∅, bijection.
- (LF) PROVED (round 123) for Construction-2.2 parents with positive homogeneous weights and condition (C):
  w(p) > w(b) for non-A p from i, b in A. (C) follows from equal A-weights (path length; |A| = 1), by induction on degree
  (d p != 0, no single-arrow word). Proof: pair each star b* with the next letter (which starts at i and has no suffix).
  Not implied by positive weights alone: kQ/(a1 c) with w(a2) = 5 gives the letter f a2^-1 of weight -3.
  Best form (scholar's review): under (C), phi(i) = c in (max w(A), min w(non-A from i)) makes EVERY kM letter positive
  (Adams-connected kM). (C) is sufficient, not necessary (kQ/(a1 c d - a2 d, a1 e), w(e) = 1). (LF) <=> some potential on all
  vertices works (Bellman-Ford); open: is (C) always reachable by reweighting? Checked: `123/theorist_shift.py`, 82/82.
  Covers loops at i and cyclic quivers. On a line the kM letter quiver is acyclic.
- On 10 LNAs (n = 4-6), all 82 vertex instances on A and A^op: (LF), (H1w), [PI], (H3) hold; H(I) = 0.
  Strict (H1) fails on 26 instances, all gate-closed under E-252 (F = Def 3.2 on A). On LNAs, a failure forces J != 0 (proof).
  Off the line it does not: the tree y:1->2, a:2->3, a2:2->4, z:3->5, rels ya, az is admitted at 2 with J = 0 and has d p = a.y.
- For no-loop parents, cancellation is E-235's anyway. The value of my route there: no omega-reduction needed, per instance.
- (H3) is the only hypothesis not implied. Observed leading words: the stripped relation's first two arrows (x_{i+1} x_{i+2}),
  or a*·p. g = path-length-like.
- Older: all 13 E-319 formal models cancel by proof; M4 and M7b need a parent automorphism (shared pivot, which also fails (H1w)).
  HH^0..2 blind on LNAs; q_S is the Cartan-level NO.

## What I tried
- 123: `rounds/123/theorist_lnaparent.py` (1 s): a Construction-2.2 builder for LNAs (Tate per cell (s, t) in order of
  width, exact Q, H = A checked); kM via 069/106; checks imported from `122/theorist_h031.py`; a vertex-general copy of
  the 106 engine (lf.run only scans range(3)!); the tree example; controls ([MM] dropped term, [LF] w(a2) = 5, planted g,
  closed ideal (a*)). The g for (H3) is found by recursion over d_H (base values for d = 0 letters), not by a box.
- 122: `rounds/122/theorist_h031.py` (13 formal models, box g-search). 118/111/109: M4/M5, M2, M1 by hand.

## Weakest points
- Bergman's diamond lemma (Adv. Math. 29, 1978) is cited from memory; E-240 is used beyond degree 2.
- The (H1w)-on-LNA and J != 0 arguments are sketches ("one Anick chain per cell and degree on a line").
- The tree parent is hand-written: d^2 = 0 is checked, H = A is not (it is the monomial Anick resolution).
- Cite E-258 (not only E-252) for F/R on rows with length 3-4 relations.
- The F/R match rests on E-252 (dimension data, length-2 relations only).

## Next
- Prove (H3) for LNA steps. Guess: the two leading shapes come from the relation's length (2 -> pivot, >= 3 -> x x top).
- A finite-dimensional loop parent (cyclic quiver, E-279) needs a non-line Construction-2.2 builder; then run (H1w)/(H3).
- Hom/End(T) for a loop model: the actual H-023 gap.
- Scholar: a published Lemma 2 for weight-graded non-complete semifree algebras (ENS math/0101145, E-320, nearly).
