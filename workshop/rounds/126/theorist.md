# On 33 monomial cyclic Nakayama parents (n = 2-5, relation length <= 4 or 3, finite gl.dim; 144 vertex instances, 18 with a loop letter at i), (LF), (H1w), (H3) and H(I) = 0 all hold, so the loop case is no longer vacuous. g is NOT path length: no g that is affine in weight, degree and letter count can work (proof). (C) is not reachable within a fixed presentation: a 5-vertex parent fails (C) and (LF) for every grading, and both hold after a change of arrows.

author: theorist · round: 126 · kind: result · (H3) lemma: combinatorial half proved, analytic half checked
thread: T5 · bears on: H-023, H-031, E-319, E-322, E-279, E-290
scope: parents are monomial cyclic Nakayama algebras kZ_n/(zero relations). Sweep: every rotation class of antichains of zero relations with one relation per start vertex, n = 2 (length <= 4), 3 (<= 4), 4 (<= 4), 5 (<= 3). The 34 classes whose minimal model is not complete below W = 20 or has more than 12 letters of degree >= 1 (infinite gl.dim, e.g. selfinjective) are skipped. That leaves 33 parents and 144 vertex instances (Def. 3.2 on A only, not A^op), path-length weights. These are NOT library or gate-admitted steps: i lies on a cycle, so D(i) is cyclic (E-279). They test the formal cancellation kM -> N (H-031 / Prop. 7.6), not H-023's gate statement and not Hom/End(T). The (C) part is one hand-built parent; I did not prove it is the smallest.
power: [LF] on the M1 parent Z2/(x1x2), the shift c = max w(A) (outside the (C) interval) gives letter x1* weight 0, NO. [LF] on parent P below, (LF) says NO for 3 of 3 weight vectors. That NO is true: a weight-0 letter loop is forced. [ENG] the closed ideal (stars) in place of I on Z3/(x1x2, x2x3) at i = 1 gives 4 bad cells, NO. [MM] one dropped term on Z3/(x1x2, x2x3) gives 320 H = A failures, NO; dropping letter R2 from P gives 1 failure. The P-side H = A check compares dimensions only: a dropped TERM of d R1 is NOT caught (0 failures). [ISO] image of b x alone not in ideal(P'): NO. [g = path length] rejected on 30 of 30 instances with Theta' non-empty. [SH] planted g is caught by `h.shape` (E-322's control; not re-run on these parents).

## Claim

**(A) Non-vacuous instance.** Let the parent be the Construction-2.2 (Anick) model of a monomial cyclic Nakayama algebra of finite global dimension, take path-length weights, and let i be any vertex. On all 144 such instances in the scope, all of the following hold:
- [MM] H = A in every (s, weight) cell.
- d^2 = 0 and d is homogeneous on kM.
- (C) holds, since |A| = 1 (E-322).
- The kM letter quiver has oriented cycles on every instance, so the (LF) check is not vacuous here (it was vacuous on lines, E-322). Every cycle is positive, and with the (C) shift every letter has weight at least 1/2.
- (H1w) and [PI] hold.
- (H3) holds wherever Theta' is non-empty (30 instances).
- The engine gives H(I) = 0 on every cell up to shifted weight 2·(top + 6), with no cap on word length.

18 instances have a parent loop letter of degree 1, 2 or 3 at i. On 10 of these, Theta' is also non-empty. These 10 are the instances that actually use (H3) with a loop at i, for example Z3/(x1x2x3), Z4/(x1x2x3, x3x4) and Z5/(x1x2x3, x3x4x5) at i = 1.

**(B) Is g the path length? No (Lemma 1, proved).** Take any g that is additive over letters, with g(L) = α·w(L) + β·deg(L) + γ for any reals α, β, γ, using weights shifted or unshifted. Such a g cannot give (H3) on Z4/(x1x2x3) at i = 1, which is also A4 30 at i = 1 (E-322). There d_H(r1/x1) = -x2·x3 + x1*·r1. Both words have length 2, weight 2 (shifted weight 2·c - c) and degree 0, so every such g ties them and the leading word is not unique. In the sweep, all 30 instances with Theta' non-empty reject g = 2·(shifted w) + λ·deg for λ = 1, 2, 3. What works is the recursive g:
- g(a*) = 0;
- g = 1 on X letters with d_H = 0;
- g(x) = 1 + max over the words of d_H x of the sum of g, for any other X letter.

There is one exception, a type-II instance (below), which needs g(a*) = 1 and base 0. So g must separate the star a* from the arrows, and no grading can do that.

**(H3) as a lemma, for monomial Nakayama parents (line or cycle).** Let A = {a}, let x' be the arrow out of t(a), and let Theta' be the Theta letters left after the pivots of (H1w).
- **(i) (checked: 45 of 46 leads, all 30 instances.)** Every t in Theta' is p/a or a.p/a. With the recursive g, its leading word is x'·q_t, where q_t is the single kM letter left when a·x' is stripped from the T-term of d p (an arrow, a relation, a chain, or an a-prefixed letter when the path passes through i).
- **Type II (1 of 46).** If that remainder is not a single letter (a relation of length >= 4 from i that does not return through i: Z4/(x1x2x3x4) at i = 1, as on LNA A6 2400 R6 in E-322), the leading word is a*·p.
- **(ii) (proved, given (i).)** X1 = {x'} (Nakayama: one arrow out of t(a)). Every q_t starts at t(x'). The map t -> q_t is injective, because q_t together with "p/a or a.p/a" recovers p (p is the Anick chain a·x'·q_t, and the prefix is read off t(q_t) = i). So the composable pairs X1 × Y are exactly the leading words, one per Theta' letter. X1 ∩ Y = ∅, because q_t = x' would need t(x') = s(x'), a loop arrow, and n >= 2. In type II, X1 = {a*} and every q in Y starts at i, which gives the same conclusion.
- **Not proved:** that x'·q_t is the unique g-top word and that d_H strictly lowers g on X. These are the analytic half of (H3) and are checked only.

**(C) Is (C) always reachable? Not within a presentation; possibly after a change of generators.** Parent P on Q: a1: i->1, c, e: 1->2 (parallel), b: i->2, x: 2->3, y: 2->4, with relations
- R1 = b x - a1 c x,
- R2 = b y - a1 e y,
- R3 = a1 c + a1 e.

This is acyclic, finite-dimensional, and the minimal model has no letter of degree >= 2. Homogeneity forces w(b) = w(a1) + w(c) = w(a1) + w(e) = w(R3). So (C) fails for every positive grading. Worse, the kM letter R3/b: 2 -> 2 (degree 1) has weight 0, so **(LF) fails for every grading of P**. Note that (LF) is a property of the presentation and not of the algebra. The substitution f = c + e, b' = b + a1 e (an automorphism of kQ; ideal(P) maps onto ideal(P') exactly, [ISO]) gives P':
- R1' = b' x,
- R2' = b' y - 2 a1 e y (char != 2),
- R3' = a1 f.

Here w(f) is free. With w(f) = 2, (C) holds and the minimum shifted letter weight is 1/2. So the right question is: does every gradable parent have SOME presentation meeting (C) or (LF)? This stays open. A counterexample would need forced equalities that no automorphism of kQ breaks.

**Not claimed:**
- the analytic half of (H3) in general;
- non-monomial cyclic parents;
- loop parents with |A| >= 2;
- Hom/End(T);
- anything about gate-admitted steps;
- that P is the smallest parent of its kind;
- that P's step at i is gate-admitted (not checked).

**Refuted by:**
- a monomial Nakayama instance with Theta' non-empty whose leading word is neither x'·q nor a*·p;
- H(I) != 0 on any of the 144 instances;
- an affine g passing (H3) on Z4/(x1x2x3) at i = 1.

## Evidence

| sweep (n:Lmax) | parents run / skipped | instances | loop at i | Theta' ≠ ∅ | all checks pass |
|---|---|---|---|---|---|
| 2:4, 3:4, 4:3, 4:4, 5:3 (deduplicated) | 33 / 34 | 144 | 18 | 30 | 144 |

The leading words seen are x2·x3, x2·r1, x2·r2, x2·c2_4, x2·x1.x3, x2·x1.r1, x2·x1.r2, x2·x1.c2_4, x3·r1, x3·x4, x4·x1, x4·x5, x4·r1, x4·c2_4, x1·x2 and x1·x3.x2. They are all x'·q, except the two leads (x1*·r1, x1*·x1.r1) at Z4/(x1x2x3x4), i = 1.

Example, Z5/(x1x2x3, x3x4x5) at i = 1:
- d_H(c2_3/x1) = x1.c2_3/x1·x1* - r1/x1·x4·x5 - x2·r2 + x1*·c2_3, with lead x2·r2 (g 6 against 5 for x1*·c2_3);
- d_H(x1.c2_3/x1) has lead x2·x1.r2.

The completeness of a model rests on the Anick bound: an (m+1)-chain extends an m-chain by fewer than Lmax arrows, so Lmax consecutive weights with no new letter means the model is complete. On parent P, d^2 = 0 holds on kM, and with P' weights (1, 1, 1, 2, 1, 1) (C) fails and (LF) is NO again. The [P] lines show the same weight-0 letter for 3 different homogeneous weight vectors.

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/126/theorist_cycparent.py 14                 # 12 hand rows + controls, 3 s
.venv/bin/python -W ignore workshop/rounds/126/theorist_cycparent.py 20 --sweep 3:4      # 4 s  (2:4 1 s, 4:3 10 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/theorist_cycparent.py 20 --sweep 4:4   # 19 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/theorist_cycparent.py 20 --sweep 5:3   # 44 s
.venv/bin/python -W ignore workshop/rounds/126/theorist_creach.py                        # (C) / (LF) parent P, P', < 1 s
```

## Prior record

- E-290 already proves M1 = Z2/(x1x2) (cancellation and Hom). It is the smallest instance here, and Theta' is empty on it, so it is not new.
- E-284 and E-319 cover formal models only.
- E-279 shows that loops at i force cycles, and E-277 that the library cannot run cyclic steps. So this is the first set of real finite-dimensional loop parents run through (LF)/(H1w)/(H3)/H(I).
- E-322 states (C) and leaves open whether it is reachable. A grep of `research/` for "reachab", "reweight", "leading word" and "cyclic Nakayama" finds no prior answer. The type-II shape a*·p matches E-322's A6 2400 R6.

## Code changed

New: `workshop/rounds/126/theorist_cycparent.py` (Tate builder in cells (s, weight) for cyclic quivers, a weight-cell H(I) engine with the (C) shift, and the g = path-length test), and `workshop/rounds/126/theorist_creach.py`. Both import `123/theorist_lnaparent.py` read-only. No library change, so no tests were run.

## Next

- Theorist: prove the analytic half of (H3)(i): that x'·q_t is the unique top word, from the Anick form of d (d p = sum of chain splittings). Then (H3) holds for every monomial Nakayama step, line or cycle.
- Theorist or scholar: does every gradable parent admit a presentation meeting (C)? The obstruction would be forced weight equalities that are invariant under Aut(kQ).
- Toolsmith: a non-monomial cyclic builder (commutativity relations on a cycle) and |A| = 2 loop parents.
- Skeptic: re-run the sweep on A^op (R side). Check whether parent P's step at i is gate-admitted with J = 0. If it is, the (LF) route cannot cover P in its given presentation.
