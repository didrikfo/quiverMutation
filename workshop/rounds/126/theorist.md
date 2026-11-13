# On all 34 monomial cyclic Nakayama parents of finite gl.dim in range (n = 2-5, relation length <= 4 or 3; 149 vertex instances, 20 with a loop letter at i), (LF), (H1w), (H3) and H(I) = 0 all hold, so the loop case is no longer vacuous. g is NOT path length: no letter-additive g that is affine in weight, degree and letter count can work (proof). (C) is not reachable within a fixed presentation: a 5-vertex parent fails (C) and (LF) for every grading, and both hold after a change of arrows.

author: theorist · round: 126 · kind: result · (H3) lemma: combinatorial half proved, analytic half checked
thread: T5 · bears on: H-023, H-031, E-319, E-322, E-279, E-290
scope: parents are monomial cyclic Nakayama algebras kZ_n/(zero relations). Sweep: every rotation class of antichains of zero relations with one relation per start vertex, n = 2 (length <= 4), 3 (<= 4), 4 (<= 4), 5 (<= 3). 68 classes in all. The 34 classes of infinite gl.dim (Kupisch-series syzygy test, `gldim_finite`) are skipped; for each, the cap-12 builder independently fails to find a complete model. The other 34 classes are all run (letter cap raised to 60; none hit it), and an independent enumeration of Kupisch series finds exactly these 34 finite-gl.dim classes (1 / 4 / 10 / 19). That gives 34 parents and 149 vertex instances (Def. 3.2 on A only, not A^op), path-length weights. These are NOT library or gate-admitted steps: i lies on a cycle, so D(i) is cyclic (E-279). They test the formal cancellation kM -> N (H-031 / Prop. 7.6), not H-023's gate statement and not Hom/End(T). The (C) part is one hand-built parent; I did not prove it is the smallest.
power: [LF] on the M1 parent Z2/(x1x2), the shift c = max w(A) (outside the (C) interval) gives letter x1* weight 0, NO. [LF] on parent P below, (LF) says NO for 3 of 3 weight vectors. That NO is true: a weight-0 letter loop is forced. [ENG] the closed ideal (stars) in place of I on Z3/(x1x2, x2x3) at i = 1 gives 4 bad cells, NO. [MM] one dropped term on Z3/(x1x2, x2x3) gives 320 H = A failures, NO; dropping letter R2 from P gives 1 failure. The P-side H = A check compares dimensions only: a dropped TERM of d R1 is NOT caught (0 failures). [ISO] image of b x alone not in ideal(P'): NO. [g = path length] rejected on 33 of 33 instances with Theta' non-empty. [GLDIM] Kupisch test on the selfinjective 2-cycle: infinite, as it must be. [KUPISCH] the sweep's run set equals the independent list of finite-gl.dim Kupisch classes for each of 2:4, 3:4, 4:4, 5:3 (missing: none); this control is the one that would have caught the round-126 hole. [SH] planted g is caught by `h.shape` (E-322's control; not re-run on these parents).

## Response to referee

Review: `theorist.review.md` (skeptic, minor revision). All four required items are done.

1. **Skip criterion; missed class; Kupisch control. Done.** The parenthetical "(infinite gl.dim)" was false: the skip was a size cap (<= 12 letters of degree >= 1). The sweep now skips a class only if the Kupisch-series syzygy test proves infinite gl.dim (`gldim_finite` in `theorist_cycparent.py`). As a cross-check it records that the cap-12 builder finds no complete model for every such class (0 exceptions). Finite-gl.dim classes run with `--cap 60`, and any class that still hit the cap would print `SKIPPED BY SIZE CAP` (none does). A new [KUPISCH] control enumerates Kupisch series independently and compares them with the run set: 1 / 4 / 10 / 19 classes for 2:4 / 3:4 / 4:4 / 5:3, with none missing. Z5/(x1x2x3, x2x3x4, x3x4x5, x5x1) is now in the sweep. All four sweeps were re-run (outputs saved, totals from `theorist_tally.py`), and they confirm the referee's numbers: 34 parents, 149 instances, 149 all-pass, 20 loop at i, 33 Theta' non-empty, 11 both, g = path length rejected on 33 of 33. On leads, my count is 52 leads: 50 are x'·q_t (32 instances) and 2 are a*·p (the one type-II instance). The referee's "51 of 52" and my earlier "45 of 46" both counted the type-II instance as one lead. Title, scope, claim (A) and evidence are updated.
2. **(B) narrowed. Done.** The claim is now "no letter-additive g affine in (w, deg, 1)", with the cases not excluded stated. I did not attempt the stronger statement.
3. **Type-II wording. Done.** The new criterion is "the remainder after stripping a·x' is not a single kM letter", and "does not return through i" is dropped. I added that the type-II branch rests on one instance.
4. **Injectivity sentence in (H3)(ii). Done (rewritten).** The a-prefix is now carried inside q_t, through |P| = a·x'·|q_t| and deg P = deg q_t + 1, and the r1/x1 vs x1.r1/x1 pair is worked. The rewrite found one more case: on Z3/(x1x2x3x1, x3x1x2), the *unprefixed* r2/x1 has q = x1.x3, so "prefix of t iff prefix of q" is false too. The argument now goes through (path, degree), with the Anick facts it relies on marked as the weakest step.

Not required, not done: the A^op check of Def. 3.2 (one line in *Next*). I did not have time this sitting to confirm Def. 3.2 on A^op against the L-side step.

## Claim

**(A) Non-vacuous instance.** Let the parent be the Construction-2.2 (Anick) model of a monomial cyclic Nakayama algebra of finite global dimension, take path-length weights, and let i be any vertex. On all 149 such instances in the scope, all of the following hold:
- [MM] H = A in every (s, weight) cell.
- d^2 = 0 and d is homogeneous on kM.
- (C) holds, since |A| = 1 (E-322).
- The kM letter quiver has oriented cycles on every instance, so the (LF) check is not vacuous here (it was vacuous on lines, E-322). Every cycle is positive, and with the (C) shift every letter has weight at least 1/2.
- (H1w) and [PI] hold.
- (H3) holds wherever Theta' is non-empty (33 instances).
- The engine gives H(I) = 0 on every cell up to shifted weight 2·(top + 6), with no cap on word length.

20 instances have a parent loop letter of degree 1, 2 or 3 at i. On 11 of these, Theta' is also non-empty. These 11 are the instances that actually use (H3) with a loop at i, for example Z3/(x1x2x3), Z4/(x1x2x3, x3x4) and Z5/(x1x2x3, x3x4x5) at i = 1, and Z5/(x1x2x3, x2x3x4, x3x4x5, x5x1) at i = 2 (loop c3_9).

**(B) Is g the path length? No (Lemma 1, proved).** Take any g that is additive over letters, with g(L) = α·w(L) + β·deg(L) + γ for any reals α, β, γ, using weights shifted or unshifted. Such a g cannot give (H3) on Z4/(x1x2x3) at i = 1, which is also A4 30 at i = 1 (E-322). There d_H(r1/x1) = -x2·x3 + x1*·r1. Both words have length 2, weight 2 (shifted weight 2·c - c) and degree 0, so every such g ties them and the leading word is not unique. In the sweep, all 33 instances with Theta' non-empty reject g = 2·(shifted w) + λ·deg for λ = 1, 2, 3. What works is the recursive g:
- g(a*) = 0;
- g = 1 on X letters with d_H = 0;
- g(x) = 1 + max over the words of d_H x of the sum of g, for any other X letter.

There is one exception, a type-II instance (below), which needs g(a*) = 1 and base 0. So g must separate the star a* from the arrows, and no letter-additive g affine in (w, deg, 1) can do that. Not shown: that no other grading can, e.g. a non-affine function of (w, deg), or a grading that gives a* a weight independent of its arrow.

**(H3) as a lemma, for monomial Nakayama parents (line or cycle).** Let A = {a}, let x' be the arrow out of t(a), and let Theta' be the Theta letters left after the pivots of (H1w).
- **(i) (checked: 50 of 52 leads, on 32 of the 33 instances.)** Every t in Theta' is P/a with P = p or P = a.p (checked on all 52). With the recursive g, its leading word is x'·q_t, where q_t is the single kM letter left when a·x' is stripped from the T-term of d P (an arrow, a relation, a chain, or an a-prefixed letter a.u when the remaining path passes through i).
- **Type II (2 of 52 leads, 1 of 33 instances).** The condition is: the remainder after stripping a·x' is not a single kM letter. Then the leading word is a*·p. Seen once: Z4/(x1x2x3x4) at i = 1 (remainder x3·x4), as on LNA A6 2400 R6 in E-322. The type-II branch of (ii) rests on this one instance.
- **(ii) (proved, given (i).)** X1 = {x'} (Nakayama: one arrow out of t(a)). Every q_t starts at t(x'). The map t -> q_t is injective. Write |q| for the parent path a kM letter stands for, with |a.u| = u followed by a. For t = P/a, |P| = a·x'·|q_t| and deg P = deg q_t + 1. So q_t fixes the path and degree of P. The a-prefix travels inside q_t, not on t's side: on Z3/(x1x2x3) at i = 1, r1/x1 has q = x3 and x1.r1/x1 has q = x1.x3 (|x1.r1| = x1x2x3x1 = a·x'·|x1.x3|), and on Z3/(x1x2x3x1, x3x1x2) at i = 1 the unprefixed r2/x1 already has q = x1.x3. Then (path, degree) determines P: two parent letters with the same path and degree are the same Anick chain, and P = p versus P = a.p' cannot share both, because then the m-chain p' would be a proper prefix of the m-chain p. *Weakest step:* the two Anick facts (chains of a monomial algebra are determined by path and degree, and no proper prefix of an m-chain is an m-chain) are cited, not reproved here, as is the identification of model letters with Anick chains. Data: no lead word repeats within any of the 33 instances. So the composable pairs X1 × Y are exactly the leading words, one per Theta' letter. X1 ∩ Y = ∅, because q_t = x' would need t(x') = s(x'), a loop arrow, and n >= 2. In type II, X1 = {a*} and every q in Y starts at i, which gives the same conclusion.
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
- H(I) != 0 on any of the 149 instances;
- an affine g passing (H3) on Z4/(x1x2x3) at i = 1.

## Evidence

| sweep (n:Lmax) | parents run / skipped | instances | loop at i | Theta' ≠ ∅ | all checks pass |
|---|---|---|---|---|---|
| 2:4 | 1 / 5 | 2 | | | |
| 3:4 | 4 / 7 | 12 | | | |
| 4:4 (contains 4:3) | 10 / 14 | 40 | | | |
| 5:3 | 19 / 7 | 95 | | | |
| total | 34 / 34 | 149 | 20 | 33 | 149 |

Every skip is infinite gl.dim (Kupisch); none is a size-cap skip. Both columns (loop at i, Theta' ≠ ∅) hold together on 11 instances. The 52 leading words, 20 distinct, are x2·x3, x2·r1, x2·r2, x2·c2_4, x2·x1.x3, x2·x1.r1, x2·x1.r2, x2·x1.c2_4, x3·r1, x3·x4, x4·x1, x4·x5, x4·r1, x4·c2_4, x4·c2_8, x4·c3_10, x1·x2 and x1·x3.x2, all x'·q with x' = x_(i+1), the arrow out of t(a) = i + 1 (50 leads), and x1*·r1, x1*·x1.r1 at Z4/(x1x2x3x4), i = 1 (2 leads). The class added in this revision, Z5/(x1x2x3, x2x3x4, x3x4x5, x5x1), has 19 letters (degrees up to 5, top weight 8), [MM] YES, every check true at all 5 vertices, ENG 0, loops at i = 2 (c3_9) and i = 5 (c3_10), Theta' non-empty at i = 1, 2, 3, and 6 leads, all x'·q.

Example, Z5/(x1x2x3, x3x4x5) at i = 1:
- d_H(c2_3/x1) = x1.c2_3/x1·x1* - r1/x1·x4·x5 - x2·r2 + x1*·c2_3, with lead x2·r2 (g 6 against 5 for x1*·c2_3);
- d_H(x1.c2_3/x1) has lead x2·x1.r2.

The completeness of a model rests on the Anick bound: an (m+1)-chain extends an m-chain by fewer than Lmax arrows, so Lmax consecutive weights with no new letter means the model is complete. On parent P, d^2 = 0 holds on kM, and with P' weights (1, 1, 1, 2, 1, 1) (C) fails and (LF) is NO again. The [P] lines show the same weight-0 letter for 3 different homogeneous weight vectors.

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/126/theorist_cycparent.py 14                 # 12 hand rows + controls, 3 s
.venv/bin/python -W ignore workshop/rounds/126/theorist_cycparent.py 20 --sweep 3:4      # 4 s  (2:4 1 s, 4:3 10 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/theorist_cycparent.py 20 --sweep 4:4   # 19 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/theorist_cycparent.py 20 --sweep 5:3   # ~105 s in parallel with the others (now includes the 14-letter class)
.venv/bin/python workshop/rounds/126/theorist_tally.py      # totals over theorist_sweep_{24,34,44,53}.out (saved outputs of the sweeps, same commands)
.venv/bin/python -W ignore workshop/rounds/126/theorist_creach.py                        # (C) / (LF) parent P, P', < 1 s
```

## Prior record

- E-290 already proves M1 = Z2/(x1x2) (cancellation and Hom). It is the smallest instance here, and Theta' is empty on it, so it is not new.
- E-284 and E-319 cover formal models only.
- E-279 shows that loops at i force cycles, and E-277 that the library cannot run cyclic steps. So this is the first set of real finite-dimensional loop parents run through (LF)/(H1w)/(H3)/H(I).
- E-322 states (C) and leaves open whether it is reachable. A grep of `research/` for "reachab", "reweight", "leading word" and "cyclic Nakayama" finds no prior answer. The type-II shape a*·p matches E-322's A6 2400 R6.

## Code changed

New: `workshop/rounds/126/theorist_tally.py` and the outputs `theorist_sweep_{24,34,44,53}.out`, `theorist_rows.out` (revision). `theorist_cycparent.py` in revision: skip by Kupisch gl.dim test, `--cap` (default 60), [KUPISCH] and tally lines. Originally: `workshop/rounds/126/theorist_cycparent.py` (Tate builder in cells (s, weight) for cyclic quivers, a weight-cell H(I) engine with the (C) shift, and the g = path-length test), and `workshop/rounds/126/theorist_creach.py`. Both import `123/theorist_lnaparent.py` read-only. No library change, so no tests were run.

## Next

- Theorist: prove the analytic half of (H3)(i): that x'·q_t is the unique top word, from the Anick form of d (d p = sum of chain splittings). Then (H3) holds for every monomial Nakayama step, line or cycle.
- Theorist or scholar: does every gradable parent admit a presentation meeting (C)? The obstruction would be forced weight equalities that are invariant under Aut(kQ).
- Toolsmith: a non-monomial cyclic builder (commutativity relations on a cycle) and |A| = 2 loop parents.
- Skeptic: re-run the sweep on A^op (R side). Per the referee, A^op of a monomial cyclic Nakayama algebra is again one with the same relation lengths, so this is redundant once Def. 3.2 on A^op at i is confirmed to be the L-side step at the image vertex (not checked this sitting). Check whether parent P's step at i is gate-admitted with J = 0. If it is, the (LF) route cannot cover P in its given presentation.
