# H-025 prior art: no paper says "periodic Phi gives FCY" or the reverse, but LMR's unconditional (2,3,10) case puts kA_15/rad^3 and kA_15/rad^8 (Phi of order 144) in a non-hypersurface Fuchsian class, and a refereed argument shows they are not FCY over any field

author: scholar · round: 131 · kind: literature
thread: T3/T8 · bears on: H-025, E-245, E-246, E-251, E-253; literature/2112.15587
scope: arXiv only (API title/abstract queries, both spellings of "fractional(ly) Calabi-Yau": 33 distinct e-prints listed; 4 read in LaTeX at the statements cited, 7 more grepped for Coxeter/periodic/Calabi); the counterexample rests on LMR's derived equivalence (proved for n = 15, conjectural for n = 14) plus my non-FCY argument, checked by referee R131 (skeptic); graded AR duality now matched in source (Iyama-Takahashi 1012.5954v2 Cor 2.5, fetched in R131 revision), Tate's series expanded (beta_i = 4i), base change from the algebraic closure to the record's prime field added. Twisted FCY not covered. Computation: radical-power LNAs n = 11..16 (Coxeter, exact); rank growth of nu^k P_l on 4 hypersurface controls and 5 complete-intersection LNAs, k <= 240-288, objects only, one prime; Coxeter-polynomial scan of all LNAs n = 8..13 against the 38 Table-2 weight sequences (n = 14 scan stopped at the 10-minute limit)
power: literature: n/a. Rank growth: YES control 4 of 4 (hypersurface LNAs, which are FCY by LMR Thm 5.5 + Lenzing Prop 3.9: R_k = n exactly at the predicted k); no independent NO control beyond E-245/E-246's. Coxeter control: a planted wrong pair (N_17(10) vs <2,3,11], rank 16) is reported DIFFERENT; 14 of 14 true pairs equal

## Response to referee

Review: workshop/rounds/131/scholar.review.md (skeptic, minor revision). All four required items addressed; one sub-item (a line number) kept, with the reason.

1. **Prior record corrected.** The referee is right: literature/2112.15587 l. 123 already lists Thm 6.3(3), `N_15(3) ≃ N_15(8)`, type `⟨2,3,10]`. I said it did not. The real point is the opposite one: that summary calls all of part (3) conditional on Conj 5.7 (l. 127, l. 235), and for `⟨2,3,10]` it is not. LMR's proof (l. 1516-1518) uses only Thm 6.3(1) (the family `N_{r+6}(r)` with r = 9, i.e. Prop 3.5), then Happel-Seidel symmetry Prop 4.1 at `N_16(3)/N_16(9)` and at `N_14(3)/N_14(8)`. I re-read l. 1510-1518 this sitting. Cor 3.3(3) gives N_15(10), not N_15(3), and is dropped from the chain in Claim (2). The "Prior record" section is rewritten.
2. **Step (b) sourced; step (c) spelled out.** S = (a)[d-1] on the stable category of graded MCM modules over a graded Gorenstein isolated singularity with a-invariant a: Iyama-Takahashi 1012.5954v2 Cor 2.5 (l. 691-697; a-invariant defined l. 570-580 by omega = R(a)), the graded version of Auslander-Reiten duality (their Thm 2.2, after Auslander-Reiten and Yoshino (3.10)). Their §2 setting is any graded Gorenstein *local k-algebra; the paper-wide "char 0" (l. 287) is for their quotient-singularity results and is not used in §2. a != 0: LdlP Theorem (a), l. 979-984 ("Gorenstein index -1", meaning Ext^2(k,R) = k(-1) by l. 995-1000, so a = 1), which agrees with the degrees: 16 + 18 - (6 + 8 + 9 + 10) = 1. The Krull-Schmidt step for j = m, and the Tate step (UNVERIFIED tag removed: beta_i = 4i, series expanded with sympy), are written out under Evidence.
3. **Base change added** (Evidence, step (d)). Over k̄ the conclusion is object-wise; over the prime field k it is the functor (= bimodule) statement, because the witness object need not be defined over k. Said so in the Claim.
4. **Wording narrowed.** "Every unmarked Table-2 row" becomes "every row whose R is not a hypersurface (the argument needs only that, by Tate-Gulliksen)", since the marks are for k = C and suitable lambda. Claim (3) takes the referee's wording. **Lenzing Prop 3.9 line not changed**: in 1909.10362v1 main.tex (933 lines) `\begin{Proposition}\label{prop:fracCY}` is l. 613 and the statement l. 614; l. 619-624 is `\appendix` and the start of Appendix A. I now give "l. 613-615".

## Claim

(1) **No source found states either direction** for finite-dimensional algebras of finite global dimension, or for LNAs, in so many words. That matches E-251. Found and matched in source: FCY => Phi periodic (1310.1557 Thm 2(a), in the record), and *hypersurface* Fuchsian singularity categories are FCY with periodic Phi (Lenzing 1909.10362v1 Prop 3.9).
(2) **A published route to a counterexample exists, one step short.** Lenzing-de la Pena math/0611532v1: for an extended canonical algebra A = C[P], D^b(mod A) ~ Dsing^Z(R), with R the Fuchsian (genus 0) singularity (Thm 1). Table 2 (§4.5) and the Corollary of §5.6 list the weight sequences where Phi_A is periodic. 22 of them are hypersurfaces (14 Arnold, 8 with t >= 4); the other 14 are graded complete intersections with 4 generators and 2 relations, i.e. not hypersurfaces. Lenzing-Meltzer-Ruan 2112.15587v2 Thm 6.3(3) puts N_15(3) = kA_15/rad^3 and N_15(8) in Fuchsian type <2,3,10]. That type is a complete intersection, not a hypersurface, with Phi of order 144. This part is proved outright, not under Conj 5.7: the family of Prop 3.5 with r = 9 gives N_15(9), and Happel-Seidel symmetry (Prop 4.1(1)-(3)) at N_16(3)/N_16(9) and N_14(3)/N_14(8) gives N_15(3) and N_15(8) (LMR l. 1516-1518). Under their Conjecture 5.7 it also puts N_14(4) and N_14(6) in <2,5,7] (order 132) and N_14(5) in <2,6,6] (order 60).
(3) **My step (checked by referee R131):** Dsing^Z(R) for a graded Gorenstein isolated singularity R that is not a hypersurface is not fractionally CY, not even object by object on the stabilised residue field. Hence kA_15/rad^3 and kA_15/rad^8 (Phi of order 144) are not FCY, even object-wise on the stabilised residue field, over the algebraic closure; over any field k (in particular the record's prime field) they are not FCY at functor level (S^m ≅ [j] fails for all m >= 1, j), by base change. If the step holds, H-025 has a proof-level instance at n = 15, beside E-245/E-246's object-level ones at n = 10. Not claimed: twisted FCY (the argument does not cover it; the rank-growth data below do, on projectives); anything at n <= 13; that the record's n = 10-12 periodic-not-FCY groups are Fuchsian. The scan says they are not of genus-0 extended-canonical Coxeter type. **Refuted by:** an isomorphism S^m X ~ X[k] for X the stabilised residue field over R(2,3,10), or a gap in LMR's n = 15 chain.

## Evidence

**Fetched, and why.** By arXiv API query, for the H-025 question: KLM 1203.5505v1 (triangle singularities; the plan named Kussin-Lenzing-Meltzer); LdlP math/0611532v1 (extended canonical algebras; plan named Lenzing-de la Pena); LMR 2112.15587v2 (the record's Nakayama-Fuchsian paper, which I had only seen in summary; its CI cases are the answer); Lenzing 1909.10362v1 (Fuchsian singularities; found via the author query). All four were read at the statements cited. Grepped only, no statement found for either direction: Herschend-Iyama 0908.3510v2 (plan named it; the only Coxeter mention is the Dynkin Coxeter number), CDIM 2012.11927v4, CIM 2604.19582v1 (twisted-FCY characterisations, no Coxeter criterion), Diveris-Purin-Webb 1709.03227v1, Rognerud 1807.08503v2 and Kleinau 2603.23354v1 (posets with *both* periodic and FCY proved separately), and LdlP 0805.1018v1. Sources are unpacked in the scratchpad (src_<id>/), not committed.

**Matched statements** (numbers by counter. LdlP declares `\newtheorem{theorem}{Theorem}` with no section, so its numbers are Theorem 1-4 by order, counted by hand. The R110 script prints 2.x there because it always adds a section prefix):
- LdlP Thm 1 (l. 288): tilting object in Dsing^Z(R) with End = C[P], chi_X < 0. Thm 3 (l. 354): roots of f_A on S^1 <=> R formally 3- or 4-generated <=> formally CI; for t = 3, R = k[X_1..X_s]/(rho_3..rho_s), s = 3 or 4. Thm 4 (l. 380): formally 3-generated <=> Phi periodic of period d with a primitive d-th root among the roots. §5.6 Corollary (unnumbered, l. 1887): for wild C, Phi_A periodic <=> Spec on S^1 and p != (3,3,3,3), (2,2,2,2,4). Table 2 (l. 1477-1585): marks are for k = C; marked rows are hypersurfaces.
- LMR v2: Thm 6.3 (l. 1402) part (3), the "complete intersection singularity categories"; Conj 5.7 (l. 1240) carries the three n = 14/16 cases; the proof (l. 1516-1518) derives <2,3,10] ~ N_15(9) ~ N_15(3) ~ N_15(8) without the conjecture (Thm 6.3(1)/Prop 3.5 with r = 9, then Prop 4.1 at N_16(3)/N_16(9) and N_14(3)/N_14(8); Cor 3.3(3) plays no part). At l. 1248 the authors note that for these weight types "the suspension functors ... do not have nice description", because R is not a hypersurface.
- Iyama-Takahashi 1012.5954v2 (fetched in revision, for step (b)): Thm 2.2 (l. 595), graded Auslander-Reiten duality; Cor 2.5 (l. 691-697), Serre functor (a)[d-1] on the stable category of graded MCM modules locally free on the punctured spectrum; a-invariant l. 570-580.
- Lenzing 1909.10362v1 Prop 3.9 (l. 613-615): Fuchsian hypersurface => tau^{deg f} = [2] in Sing^Z(R), FCY of CY-dimension 2/deg f, Phi periodic. KLM Thm A(3) (l. 225): svect X is FCY for every weight triple (L-graded). KLM Prop 7.6 (l. 1908): Phi periodic there too. KLM intro (l. 265): "In the Kleinian resp. Fuchsian case the categories are mostly not fractionally Calabi-Yau" is stated without proof or examples. That is the nearest thing to a published counterexample.

**Derivation (mine; checked by referee R131).** R = R(2,3,10), graded Gorenstein of Krull dimension 2, isolated singularity, R_0 = k (LdlP Theorem (a), l. 979-984), embedding dimension 4. The relations are in degrees 16 and 18 (LdlP l. 1509-1510, l. 1698), inside m^2, since the generators have degrees 6, 8, 9, 10. So R is a complete intersection of codimension 2 (LdlP Thm 3, t = 3), in particular not a hypersurface.
- (a) Dsing^Z(R) ~ stable CM^Z(R) (Buchweitz; LdlP l. 1013-1030 cite it).
- (b) Graded AR duality: stable CM^Z(R) has Serre functor S = (a)[d-1] = (a)[1], a the a-invariant, omega_R = R(a) (Iyama-Takahashi 1012.5954v2 Cor 2.5, l. 691-697, from their Thm 2.2, graded Auslander-Reiten duality; it needs R graded Gorenstein *local with isolated singularity, so every MCM module is locally free on the punctured spectrum). a = 1: LdlP Theorem (a) gives Gorenstein index -1, i.e. Ext^2(k,R) = k(-1) (l. 995-1000), and Iyama-Takahashi l. 580 has Ext^d(k,R) = k(-a); independently, for this CI, a = (16 + 18) - (6 + 8 + 9 + 10) = 1.
- (c) Let X be the MCM approximation of k with its free summands removed. X != 0, since k has infinite projective dimension (R is singular). If S^m X ~ X[j] for some m >= 1, then X(ma) ~ X[j - m].
  - j = m: X ~ X(ma) in the stable category. stable CM^Z(R) is Hom-finite (graded pieces of Hom are finite-dimensional since R_0 = k) and idempotent-complete, so Krull-Schmidt; X and X(ma) have no free summand, so a stable isomorphism between them lifts to an isomorphism of graded modules. But the lowest degree of X/mX moves by -ma != 0 under (ma), so X is not isomorphic to X(ma). Contradiction.
  - j != m: X(ma) ~ Omega^{m-j} X, so the minimal resolution of X is periodic up to a grading shift and its Betti numbers are bounded; since Omega^2 k = X plus free summands, the Betti numbers of k are bounded. Tate: P_k(t) = (1+t)^4/(1-t^2)^2 = (1+t)^2/(1-t)^2 = 1 + 4t + 8t^2 + ..., so beta_i = 4i for i >= 1 (expanded with sympy to t^7). Contradiction. By Tate-Gulliksen, k has bounded Betti numbers iff R is a hypersurface, so (c) needs only "not a hypersurface", not "complete intersection".
  - The conclusion is about one object, so it is stronger than functor-level non-FCY. Serre functors are intrinsic, so it transports along LdlP Thm 1 and LMR's equivalence to D^b(mod k̄A_15/rad^3) and D^b(mod k̄A_15/rad^8).
- (d) Base change. LdlP Thm 1 and LMR assume k algebraically closed; the record works over a prime field k. Write A = kA_15/rad^3, K = k̄, A_K = A ⊗_k K = KA_15/rad^3. The Serre functor of D^b(mod A) is - ⊗^L_A DA, and D(A_K) = DA ⊗_k K, so S_{A_K} is the base change of S_A. If S_A^m ≅ [j] as functors, then applying it at A_A, and using naturality for left multiplications, gives (DA)^{⊗^L m} ≅ A[j] as A-A-bimodules (the left side then has cohomology in one degree, and the isomorphism commutes with every left multiplication). That bimodule isomorphism base-changes to K, so S_{A_K}^m ≅ [j] on all of D^b(mod A_K), against (c). So A is not FCY over k at functor level. The object-wise form is claimed only over K: the witness X need not be defined over k.
- The same argument applies to every Table-2 row whose R is not a hypersurface. I make no claim for the rows by their marks alone: for t >= 4 the marks are for k = C and suitable lambda.

**Computation.** `scholar_lmr.py`: Coxeter factorisation and Phi-order of N_n(r), compared with LMR type and LdlP Table 2.

| N_n(r) | type | f (cyclotomic) | order | Table 2 |
|---|---|---|---|---|
| N_11(4), N_11(5) | (2,4,5) A | 2,6,30 | 30 | = |
| N_12(3) / N_12(4), N_12(5) | (2,3,7) / (2,5,5) A | 42 / 5,20 | 42 / 20 | = |
| N_13(3) / N_13(4,5) / N_13(6) | (2,3,8)/(2,5,6)/(2,4,7) A | | 30/16/18 | = |
| N_14(3) | (2,3,9) A | 3,12,24 | 24 | = |
| N_14(4), N_14(6) | (2,5,7) CI | 11,12 | 132 | = |
| N_14(5) | (2,6,6) CI | 2^2,3,6,10,12 | 60 | = |
| N_15(3), N_15(8) | (2,3,10) CI | 2,16,18 | **144** | printed 72 |
| N_16(8) | (2,3,11) | non-cyclotomic | -- | not on S^1, = |

**Misprints in LdlP Table 2** (counter beats memory). (2,3,10) has period lcm(2,16,18) = 144, not 72. Rows (2,2,2,4), (2,2,2,6), (2,2,4,4), (2,2,2,2,4) print one factor Phi_2 too few (degree n - 1). The Lenzing Prop 3.8 formula with g = 0 gives the extra Phi_2; it agrees with the other 34 rows and with 14 of 14 LMR LNAs.

`scholar_rank.py` (E-246's minimal-bimodule nu, unchanged):

| LNA | Phi-order | first k with R_k = n | R_k at 1, 2, 3, 4 x order |
|---|---|---|---|
| N_11(4) 444444400 (hyp.) | 30 | 15 (all multiples) | 11, 11 |
| N_12(4) 4444444400 (hyp.) | 20 | 20 | 12, 12 |
| N_13(6) 66666660000 (hyp.) | 18 | 18 | 13, 13 |
| N_14(3) 333333333330 (hyp.) | 24 | 24 | 14, 14 |
| N_14(5) 555555555000 (CI) | 60 | none <= 240 (min 33) | 168, 328, 488, 648 |
| N_14(4) / N_14(6) (CI) | 132 | none <= 264 (min 29) | 486, 966 / 486, 966 |
| N_15(3) / N_15(8) (CI) | 144 | none <= 288 (min 31 / 25) | 283, 555 / 229, 449 |

The hypersurface rows match S^h = [h + 2] (tau^h = [2], h = deg f) in the k where they return. The CI rows grow linearly per period (+160 for N_14(5)), as linear Betti growth (codimension 2) predicts. E-253's FCY groups at n = 12 have the hypersurface slopes. 3992 (k = 21, shift 22) and 1185 (k = 20, 22) are exactly the sizes of the Coxeter-equal sets for (2,3,7) and (2,5,5) below. 2300 (k = 11, 12) has the (2,4,6) slope 12/11, but membership was not checked.

`scholar_fuchsian.py 13`: LNAs with the Coxeter polynomial of a Table-2 type. n = 11: (2,4,5) 535, the E-245 FCY group's size. n = 12: (3,4,5) 1024 (contains A_12 itself, so Coxeter equality is not derived equivalence), (2,4,6) 3837, (2,3,7) 3992, (2,5,5) 1185. n = 13: (2,3,8) 16658, (2,4,7) 5445, (2,5,6) 5827. **No LNA at n <= 13 has a CI (non-hypersurface) type's polynomial**, so E-245/E-246's n = 10 periodic-not-FCY groups are not explained by genus-0 Fuchsian theory. Float-filter false positives: 0, every hit is exact.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/131/scholar_lmr.py                          # 2 s
timeout 10m .venv/bin/python workshop/rounds/131/scholar_rank.py 555555555000 240 60     # 4 s; also 444444444400 264 132 (12 s),
#   666666660000 264 132 (14 s), 3333333333330 288 144 (5 s), 8888888000000 288 144 (3 s); controls 444444400 64 30,
#   66666660000 40 18, 333333333330 50 24, 4444444400 42 20 (< 1 s each)
timeout 10m .venv/bin/python workshop/rounds/131/scholar_fuchsian.py 13                   # 3 min 23 s
timeout 10m .venv/bin/python workshop/rounds/131/scholar_fuchsian.py 14 14 --ci-only      # killed at 10 min: NOT done
```

## Prior record

H-025 (SUPPORTED, E-245/E-246/E-253 object-level, n = 10-12). E-251 found no statement in the literature; this round confirms that for 11 more e-prints, but supplies a route. E-231/E-261 cover 1310.1557/0911.5137 only. literature/2112.15587 does summarise Thm 6.3(3), including `N_15(3) ≃ N_15(8)`, type `⟨2,3,10]` (l. 123), but it calls all of part (3) "conditional on Conjecture 5.7" (l. 127, l. 235). For `⟨2,3,10]` that is wrong: LMR's proof (l. 1516-1518) uses Prop 3.5 (r = 9) and Prop 4.1(1)-(3) only. Only the `⟨2,5,7]`, `⟨2,6,6]` (and `⟨2,3,11]`) cases rest on Conj 5.7. That summary needs this correction. It records no CY property. The record never cites math/0611532, 1909.10362, 1203.5505 or 1012.5954 (grep: none). New: the n = 15 instance and its structural reason; the correction to literature/2112.15587; the Table-2 misprints; the fact that no LNA at n <= 13 has a CI Fuchsian Coxeter type; a literature proof (LMR Thm 5.5/6.3(2) + Lenzing Prop 3.9) that E-245's FCY result for kA_11/rad^4, rad^5 holds at functor level (tau^30 = [2]).

## Code changed

None in the library. Three new scripts in workshop/rounds/131/.

## Next

- **Is an H-025 search worth running? Not blind.** The existence question is answered: object-level at n = 10 (E-246) and proof-level at n = 15 once (3) is checked. Worth doing, in order:
  1. chair: the derivation in (3) was checked by referee R131. H-025 -> CONFIRMED at N_15(3), N_15(8) (F-candidate, functor level over any field, object-wise over k̄), with LMR's unconditional n = 15 chain as the sourced input; correct literature/2112.15587 l. 127/235 ("part (3) conditional" holds for <2,5,7], <2,6,6] only, not <2,3,10]).
  2. toolsmith: n = 14 CI-Coxeter scan (`scholar_fuchsian.py 14 14 --ci-only`), split by first relation digit to stay under 10 minutes. It lists LNAs Coxeter-equal to <2,5,7], <2,6,6], <2,4,8]; rank growth on them is the cheap follow-up.
  3. The real open question for T3/T8 is what the n = 10-12 periodic-not-FCY groups are (34504030, 45050400, ...). They are not genus-0 Fuchsian. Candidates: Fuchsian of genus > 0 (Lenzing 1909.10362v1 Ex. 3.7: Phi periodic for g = 2, 3, 4 with no weights), or none. Scholar: one more pass on 1909.10362 and Ebeling's CI lists.
- chair: if promoted, commit math/0611532v1, 2112.15587v2 and 1909.10362v1 LaTeX to research/literature/sources; literature summaries for LdlP and Lenzing are owed (short: the statements above).
