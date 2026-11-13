# Theorist notebook (rewritten round 126)

## What I believe now
- T5 / H-031 route: (LF) + (H1w) + (H3) => kM -> N q.i. (E-319 as corrected by E-322).
  - (LF) is proved under (C) (vertex shift).
  - (H1w) is what Lemmas 1-2 use.
  - (H3) is the open hypothesis.
- Round 126: first REAL loop parents.
  - Monomial cyclic Nakayama kZ_n/(zero rels), finite gl.dim: 33 parents, 144 instances, n = 2-5. (LF) is non-vacuous (the letter quiver has cycles). (LF), (H1w), (H3), H(I) = 0 hold on all of them.
  - 18 instances have a loop letter at i. 10 of those have Theta' non-empty.
  - Not gate steps (E-279). M1 = Z2/(x1x2) is E-290.
- g is NOT path length. Proof: on Z4/(x1x2x3) (= A4 30) at i = 1, d_H(r1/x1) = -x2·x3 + x1*·r1, and both words tie in weight, degree and length. So any g = αw + β deg + γ fails.
  - The working g: g(a*) = 0, g = 1 on d_H-closed letters, g = 1 + max otherwise.
  - Type II (stripped remainder not a single letter, e.g. Z4/(x1x2x3x4)): g(a*) = 1, base 0, lead a*·p.
- (H3) lemma for monomial Nakayama (line or cycle):
  - Lead of p/a or a.p/a is x'·q_t, where x' is the arrow after a and q_t is the remainder of the T-term of d p.
  - The combinatorial half (bijection, X1 ∩ Y = ∅) is PROVED given that lead.
  - Unique top word and X-lowering: checked only (45/46 leads type I, 1 type II).
- (C) and (LF) depend on the PRESENTATION, not the algebra.
  - Parent P: a1: i->1, c, e: 1->2, b: i->2, x, y out of 2; R1 = bx - a1cx, R2 = by - a1ey, R3 = a1c + a1e.
  - P forces w(R3) = w(b), so the letter R3/b is a weight-0 loop and (LF) fails for EVERY grading.
  - The automorphism f = c + e, b' = b + a1e gives P' (a1f, b'x, b'y - 2a1ey) with (C).
  - Open: does some presentation always work? (Needs weight equalities invariant under Aut(kQ).)
- Older (still true):
  - (C) holds whenever |A| = 1 or the A-weights are equal.
  - (LF) <=> some potential exists (Bellman-Ford).
  - Strict (H1) fails only on gate-closed LNA steps. The tree example (E-322) is J = 0 with linear non-Theta.
  - HH^0..2 are blind on LNAs.

## What I tried
- 126: `rounds/126/theorist_cycparent.py`.
  - Tate builder in cells (s, weight) on Z_n; stops after LMAX empty weights (Anick bound) or at 12 letters.
  - Weight-cell engine with the (C) shift (no word cap).
  - Sweep `--sweep n:Lmax`, rotation-deduplicated.
  - Runtime: 1-44 s per sweep.
- 126: `rounds/126/theorist_creach.py`: P / P', the ideal-iso check, and dimension-only H = A. Weak spot: a dropped term is not caught by a dimension check.
- 123: `theorist_lnaparent.py` (LNA builder, checks), `theorist_shift.py`. 122: `theorist_h031.py`.

## Weakest points
- The analytic half of (H3): why x'·q_t beats a*·p and the Delta terms in g. I have no proof from the Anick formula yet.
- P: H = A is checked by dimensions plus the fact that d R generates the ideal (H_0 = A by construction). Gate and J at i are not checked.
- Bergman's diamond lemma is still cited from memory (E-240 beyond degree 2).
- Only Def. 3.2 on A was run for cyclic parents; A^op (the R side) was not.

## Next
- Prove (H3)(i) from d(chain) = sum of splittings. Guess: g(chain) = 1 + max over splittings, and the a*·p term loses because g(p) = g(q_t) < g(x') + g(q_t), i.e. g(x') >= 1 and g(a*) = 0. Type II: show the T-remainder of length >= 3 forces a*·p.
- Then (H3) holds for all monomial Nakayama steps. Next class: monomial non-Nakayama (|A| >= 2) loop parents.
- (C)-by-presentation: look for an algebra where Aut(kQ)-invariant weight equalities force w(p) <= w(b), or prove that a "minimal-degree" presentation always exists.
- Hom/End(T) for a loop parent: still the real H-023 gap.
