# H-031 as three lemmas, on 13 formal dg loop models only: with a positive-cycle condition (LF), the substitution and trivial-part steps are proved with no completion, no direct-sum splitting and no product pivot; Lemma A then applies to all ten E-284 models plus M4', M5', M9, so M3, M3b, M7, M7b, M8A, M8B and the new M9 now cancel by proof

author: theorist · round: 122 · kind: proof
thread: T5 · bears on: H-031, H-023, E-284, E-290, E-294, E-309, E-315, E-240
scope: formal dg models only: the ten E-284 loop parents (M1-M5, M3b, M7, M7b, M8A, M8B) plus three new ones (M4' = M1 + a parallel a2; M5' = M1 + spectators m, e; M9: d f = ab + mec, c != b), 13 in all, printed Def. 3.2 differentials with signs `last = -1` (E-202). No Construction-2.2 parent of a finite-dimensional algebra, no Hom/End(T), so H-023 is not touched. The g-search is a box search (values 0..5 or 0..3): "no g found" would mean "none in the box", but every model in the list got one. On M9 the claim is only that E-309's product-pivot recipe does not apply; whether some other triangular substitution brings M9 to H-031's shape is not decided.
power: proof (n/a) for Lemmas 1-3. Each machine check says NO on a control where NO is right: [LF] toy quiver with a weight-0 cycle: NO; [PI] planted sign in one pivot solution: 2 dg-map failures, NO; [PSI] one dropped term in the letter map (M5, M7b): NO; [SH] planted g (largest value set to 0): NO on 10 of 10 models, but this tests only the clause "d_H strictly lowers g on X"; [SH] near-miss ideals with H(I) != 0, where (H1) and [PI] pass: I + (a*) of M1 (5 cells), I + (f) of M2 (engine W = 9: 2 bad cells), I + (a*) of M2 (12 bad cells), I + (a*) of M7 (7 bad cells): no Lemma A shape in the box, NO on all four. These are the only controls that test the shape clauses together. [H1] says NO on M4 and M7b, where NO is not the truth (both cancel); that NO means "needs a linear change first", so it is a limit of the check, not power.

## Response to referee

1. *M3 (LF), title count, g column.* Done. M3's minimum cycle weight is 2 (printed `{0: 2, 1: 2}`), not 3. The title now says 13 models (10 + 3). The g column now lists every nonzero value the script prints, Theta' letters included (e.g. M5: a.f/m 1, m.f/m 1, f/m 3; M9: f/m 2, a.f/m 1, m.f/m 1).
2. *(LF), (H1), [PI] printed for every named model.* Done. `part_models` now includes M2, M3, M8A, M8B and M7b. Printed: (LF) holds on all five (minimums 3, 2, 3, 3, 3); (H1) and [PI] pass on M2, M3, M8A and M8B; (H1) says NO on M7b (shared pivot), as stated. These agree with the referee's run.
3. *Pivot coefficient.* Done, in the lemma's direction. [H1] now accepts any c != 0, and `projection` already divides by c (`Fraction(-k, c)`). No verdict changes: every pivot coefficient in these 13 models is ±1, so the wider test is not exercised by a non-unit c here.
4. *Near-misses with H(I) != 0.* Done. `part_nearmiss` now adds M2 + (f), M2 + (a*) and M7 + (a*). On each one, (H1) and [PI] pass, [SH] finds no shape (256 g tried, NO), and the engine at W = 9 gives 2, 12 and 7 bad cells (sumH 2, 14, 19). The *Power* line now says that the planted-g control tests only the lowering clause.
5. *M9 against H-031's wording.* Narrowed as the referee proposed (Claim and *Scope*). M9 lies outside E-309's product-pivot recipe, because no a -> a + path absorbs mec. Whether some other triangular substitution brings M9 to H-031's shape is not decided. M9 cancels through (H3) with X1 = {e}. I have not shown that no such substitution exists: that needs a search over all triangular automorphisms, and I have neither an argument nor a bounded search for it in this sitting.
6. *E-240 beyond degree 2.* Cited under weakest step (2). This is the referee's run, not mine, and I have not re-run it.

The referee's extra models (M9w1, M9w2, M9coef2, M59, M9ee, M5w, M7w; no counterexample) and the smaller-ideal YES cases are the referee's evidence. They are not added to my script, and I do not claim them as mine.

## Claim

Setting: kM = T_R(V), R = k^n, on finitely many letters, each with a source, target, degree and weight; d has weight 0 and d^2 = 0; Theta = letters with a suffix (phi b^-1); I = (Theta, d Theta); N = kM/I. **(LF)**: every oriented cycle in the letter quiver of kM has positive weight. Then each (src, tgt, weight) cell is finite and word length in it is bounded.
**Lemma 1 (substitution; proved, elementary).** Assume (LF). Take an R-algebra endomorphism that fixes the idempotents and sends each letter x to x plus words of length >= 2 in x's cell (or to an invertible linear combination within a cell). It is an automorphism: it is unipotent for the finite length filtration of each finite cell. *Pivot form:* suppose **(H1)**: only Theta letters have a linear term in d; no linear term is a Theta letter; and the pivot letters x_theta (d theta = c x_theta + ..., c != 0) are distinct and occur in no other theta's linear part. Then replacing each x_theta by d theta gives free generators H ⊔ Theta_piv ⊔ d Theta_piv. Unlike DWZ Lemma 4.7, no limit is needed: (LF) makes every cell its own completion.
**Lemma 2 (trivial-part removal, non-completed; proved).** Assume (LF) and (H1). Let E = span(Theta_piv, d Theta_piv). Then (E) is a dg ideal, pi: kM -> T(H) = kM/(E) is a dg map, and pi is a quasi-isomorphism. Proof: filter by the number of E-letters (in the new generators). d does not lower that number. gr = (T(H), d_H) * (T(E), d_E) with E contractible, so gr^{>=1} is a direct sum of tensors with a contractible factor over semisimple R, hence acyclic. By (LF) the filtration is finite in each cell, so (E) is acyclic. B * T(E) does **not** have to be a dg splitting: d of H-letters may have E-terms. This is the step KY Lemma 2.10 does by Mittag-Leffler. (E) ⊆ I, and pi(I) = (Theta', d_H Theta') with Theta' = Theta \ Theta_piv.
**Lemma 3 (Lemma A on the reduced part; proved given E-240 and Bergman's diamond lemma).** Let X = H \ Theta'. Suppose **(H3)**: there is g: H -> Z>=0 such that d_H strictly lowers g on every x in X; for each theta in Theta' the g-top part of d_H theta is one monomial ±x_theta y_theta with x_theta, y_theta in X; X1 = {x_theta} and Y = {y_theta} are disjoint; and theta -> (x_theta, y_theta) is a bijection onto the composable pairs of X1 x Y. Then T(H) -> T(H)/(Theta', d_H Theta') is a q.i. Proof: the g-filtration is exhaustive and bounded below, so a map that is a q.i. on gr is a q.i.; no cell finiteness is needed. gr T(H) has d0 theta = ±xy and d0 = 0 on X, which is E-240's setting (spectators and loops allowed). gr J = (Theta', x_theta y_theta): modulo (Theta') the relations are x y - (lower g), and their leading words x y have no overlaps because X1 ∩ Y = ∅, so the diamond lemma says the normal words form a basis. The lower-order terms need **not** lie in (Theta'). The round-118 proof assumed they did; M5 and M9 violate that.
**Theorem.** (LF) + (H1) + (H3) imply that kM -> N is a q.i. (pi, then Lemma 3; N = T(H)/pi(I)).
**Tests.** (LF), (H1) and the pi-checks hold, and an (H3)-g exists, on M1, M2, M3, M3b, M5, M7, M8A, M8B, M4', M5', M9. All of these are printed by the script. (LF) is also printed for M4 and M7b. For M4 and M7b, (H1) holds only after a parent automorphism: M4 = Psi(M4') by a1 -> a1 - a2, M7b = Psi(M7) by h -> h + f. Each Psi is checked as a dg isomorphism with Psi(I') = I. **M5 (E-309)** passes directly, with g(m*) = 1, g(fbar) = 2, g(t_m) = 3. So E-309's product-pivot step was never needed. It is the image of the parent automorphism a -> a + me (M5 = Psi(M5'), checked), which E-309 did not use. **Loop case: M3b** (degree-1 f with d f = ab, degree-2 g with d g = f ab - ab f; E-284's two-loop model; H-031's hypotheses do not cover it). Pivots f/a -> b, a.f/a -> a.b. Theta' = {g/a, a.g/a}. g = (fbar 1, a.g 2, gbar 3) gives leading words a*·gbar, a*·a.g, with fbar and a.f as spectators. **New proofs:** M3, M3b, M7, M7b, M8A, M8B and M9 now cancel for all weights (previously E-284 computation to weight 9-12). **M9 is outside E-309's product-pivot recipe** (no a -> a + path absorbs mec, since c is not b). Whether some other triangular substitution brings M9 to H-031's shape is not decided, so M9 is *not* claimed as a counterexample to H-031's wording. M9 cancels through (H3) with X1 = {e}. Not claimed: anything about a library step (E-279), Hom/End(T), or that (H3) holds for every loop model. **Refuted by:** a model with (LF) + (H1) + (H3) and H(I) != 0 (an error in a lemma), or an error in E-240.

## Evidence

`theorist_h031.py` (about 1 s). [PI] verifies that pi d x = d_H pi x on **every** letter, and that d_H^2 = 0. [H1] accepts a pivot coefficient c != 0 (the projection divides by c). [SH] searches g and checks (H3) mechanically: X letters strictly lowered; unique quadratic X-monomial top; disjointness; bijection with composable pairs. [PSI] checks d Psi = Psi d on all letters, Psi(Theta) ⊆ (Theta) and Psi^-1(Theta) ⊆ (Theta), with Psi^-1 found triangularly. Injectivity follows from equal letter counts per cell and (LF). As a cross-check, the engine (E-284's `lf.run`) gives H(I) = 0 for M5' (72 cells, dim I 4 452) and for M9 (72 cells, dim I 7 516) to weight 12.

| model | (LF) min cycle wt | pivots | Theta' | g found (every nonzero value printed, Theta' included) | leading words (X1 = ...) |
|---|---|---|---|---|---|
| M1 | 3 | f/a->b, a.f/a->a.b | none | - | Lemma 3 vacuous (E-290's case) |
| M2 | 3 | none | f/a, a.f/a | b 1, f 1, f/a 1 | a*·f, a*·a.f (= E-294's g) |
| M4 | 3 | (H1) NO: shared pivot | | via M4' and Psi | a1 -> a1 - a2 checked |
| M4' | 3 | 3 | 3 | f 1, f/a2 1 | a2*·{f, a1.f, a2.f} |
| M5 | 4 | f/a->b, a.f/a, m.f/a | f/m, a.f/m, m.f/m | m* 1, f 2, f/m 3, a.f/m 1, m.f/m 1 | m*·{f, a.f, m.f}; e·a*·f is lower |
| M5' | 4 | same | same | f 1, f/m 1 | m*·{f, a.f, m.f}; M5 = Psi(M5') |
| M9 (new) | 4 | same 3 | 3 | f 1, c 2, a.c 1, m.c 1, f/m 2, a.f/m 1, m.f/m 1 | e·{c, a.c, m.c}; m*·f etc. lower |
| M3 | 2 | none | 4 | b 1, f 1, a.g 2, g 3, f/a 1, g/a 3, a.g/a 2 | a*·{f, a.f, g, a.g} |
| M3b | 3 | 2 | g/a, a.g/a | f 1, a.g 2, g 3, g/a 3, a.g/a 2 | a*·{g, a.g} |
| M7 | 3 | 2 | h/a, a.h/a | f 1, h 1, h/a 1 | a*·{h, a.h} |
| M7b | 3 | (H1) NO: shared pivot | | via M7 and Psi | h -> h + f checked |
| M8A / M8B | 3 | 2 | 4 | f 1, a.g 2, g 3, a.k 3/4, k 4/5, g/a 3, a.g/a 2, k/a 4/5, a.k/a 3/4 | a*·{g, a.g, k, a.k} |

**Weakest steps.** (1) Lemma 3 uses Bergman's diamond lemma (Bergman 1978, Adv. Math. 29). This is from memory: theorem number not matched, and its path-algebra form over R is assumed. (2) E-240 is used with spectators, loops and degrees up to 4. Its proof is degree-free and allows loops; its own machine check covered degrees <= 2 only. The referee ran E-240's `check` (`086/theorist_lemmaA.py`) on 120 random quivers (n <= 3, loops and parallels, degrees 0..4, weight <= 5): 5 265 cells, 0 mismatches, plus a 50-cell vertex-i loop shape, 0 mismatches. This supports "degree-free" but does not prove it. (3) The signs are the printed `last = -1` throughout. (4) Lemma 2 needs (LF), which is checked per model, not proved for Construction-2.2 parents.

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/122/theorist_h031.py              # all parts and controls, incl. 3 engine near-misses at W = 9, about 1 s
.venv/bin/python -W ignore workshop/rounds/122/theorist_h031.py --engine 12  # adds E-284 engine H(I) for M5', M9, about 2 s
```
I read the assignment's "E-203/E-240" as E-284's loop list: E-203 is an n = 9 T10 control with no loop, and E-240 is Lemma A itself. M3b is the loop case I worked by hand, as the next agenda item (H-023/H-031 M3/M3b).

## Prior record

E-284 (computation to weight 12/9), E-290 (M1), E-294 (M2, same g as found here), E-309 (M4, M5 with product pivot), E-315 (no published lemma; mechanism only), E-240 (Lemma A). A grep of `research/` for automorphism/naturality, locally finite, minimal model, diamond/Bergman found nothing relevant. New in this submission: (LF) as the substitute for completion; Lemma 2 without a dg splitting; Lemma 3 with lower terms outside (Theta'); the parent-automorphism reading of the product pivot; proofs for M3, M3b, M7, M7b, M8A, M8B, M9. Not new: the three-step decomposition (H-031's text), and M1, M2, M4, M5 (already proved; reproved here more simply). This is close to the standard minimal-model splitting of semifree dg algebras (the scholar's Van den Bergh / Kalck-Yang leads). If one of those states Lemmas 1-2 over R with a weight grading, cite it instead.

## Code changed

New script `workshop/rounds/122/theorist_h031.py`. It imports the round-106 scripts read-only. Revised after review: [H1] uses c != 0; `part_models` covers all 13 models; `part_nearmiss` adds three near-misses with engine sumH. No library change, no tests run.

## Next

- Keeper: H-031 should be restated as "(LF) + (H1) => trivial part cancels; cancellation holds iff/when (H3)". The open question becomes whether (H3) holds for every loop model. Falsifier: a model with (LF) + (H1) and no (H3)-g that has H(I) != 0, or one with no g where H(I) = 0 (then (H3) is too strong).
- Skeptic: second implementation of [PI]/[SH] on random loop parents (`069/scholar_oppermann_sign.random_instance` plus weights). Count how often (LF), (H1) and (H3) hold, and run the engine where they fail. Also re-derive Lemma 3's gr J argument.
- Scholar: match Bergman's diamond-lemma theorem number. Check whether Van den Bergh 1008.0599 or Kalck-Yang state Lemma 2 for weight-graded non-complete semifree algebras.
- Theorist: (LF) for Construction-2.2 parents (is it implied by positive parent weights?); (H1) when other relations pass through b; Hom/End(T) for a finite-dimensional loop model.
