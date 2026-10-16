# Oppermann v1 Sect. 8: Lemma sign is right as stated (l. 1058), the proof's last line (l. 1105) is mis-simplified, d(alpha)=0 is automatic, and a d^2 test passes

## Response to referee

1. Citation: done. "A = degree-0 arrows out of i, none a loop" is at l. 321 (l. 308 was the Remark on R End); fixed in Claim 2, Claim 3(b) and Next. l. 262 (non-negative grading) unchanged.
2. Seed and skeleton: done, in the Scope line (sign script seed 0, lemma script seed 1, skeleton listed). Not done: the split by |A| = 1, 2, 3 (would need a new tally; one line, skipped as not changing the claim). Referee's 1000- and 3000-trial runs (0 failures printed, 4019 and 12189 control) are noted; same stream prefix, so not independent.
3. Shared convention: stated in Scope. Hand check of `M.deg` against the .tex (l. 332, 342, 351, 358): ||phi|| = |phi| - chi(s phi) + chi(t phi); ||alpha phi|| = ||phi|| - 1; ||phi alpha^-1|| = ||phi|| + 1; ||alpha phi beta^-1|| = ||phi||; ||alpha^*|| = 0. `M.deg` is `dp - (sp==i) + (tp==i)`, then -1 if a, +1 if b, 0 for alpha^*: identical on all five cases. The test therefore uses the paper's own convention; remaining shared-reading risk is only my transcription of Def. 3.2 itself.
Claim, title and result are unchanged by these answers.


author: scholar · round: 069 · kind: result
thread: T5 · bears on: H-015 (J = 0 premise; no direct bearing), E-189, E-196
scope: arXiv 1504.02617v1 LaTeX only (`research/literature/sources/1504.02617v1-oppermann-quivers-silting-mutation.tex`); my own implementation of Def. 3.2 (l. 383-391) tested on 200 random instances (`random.Random(0)`, the default seed; fixed skeleton of `random_instance`: vertices 0=i,1,2,3; A = {0->1} plus 0->2 and a parallel 0->1 each with prob. 1/2; 1-2 degree-0 arrows from {1,2,3} into i; fixed degree-0 path 1->3->2; 1-3 degree-1 loops at i; 1-3 further degree-1 arrows; d of degree-1 arrows = up to 3 random degree-0 paths of length <= 4; up to two degree-2 arrows; 7 082 M-letters; lemma script seed 1). Both scripts share one M-degree function (`M.deg`), so they test the signs given that reading of ||.||, not the reading itself. Not covered: arrows of degree >= 3, journal version (UNVERIFIED), the quasi-isomorphism proof (Sect. 7).

## Claim

1. Sign (E-196, first item). The Lemma statement (l. 1055-1059), second formula, is correct as printed (modulo a typo: `\ta{\varphi}` should be `\ta{p}`, i.e. t(p)):
   d(Δp α^-1) = d(Δp) α^-1 + (-1)^|p| Δp α^* - (-1)^{c_i(t p)} Δ(dp/α), with |p| meaning the M-degree ||p||.
   The last line of its proof (l. 1105) reads "+(-1)^{||p_2||} dec(dp)/alpha" and does not match; the underbrace at l. 1104 "p_2 (d p_1)/alpha = (-1)^{||p_2||} dp/alpha" is wrong as printed. Leibniz for d on kQ gives (dp)/α = (-1)^{|p_2|} p_2 (dp_1/α) (the term dp_2 p_1 has no path starting with α because p_1 does not start in A), with |p_2| the Q-degree, and ||p_2|| = |p_2| - 1 + c_i(t p_2) since s(p_2) = i. So the two signs differ by -(-1)^{c_i(t p)}, and the correct last term is (-1)^{||p_2||}(-1)^{|p_2|} Δ(dp/α) = -(-1)^{c_i(t p)} Δ(dp/α), which is the statement. Also the inductive input at l. 1100 (written "+") is right because in that case t(p_1) = i, so -(-1)^{c_i(t p_1)} = +1. Net: a typo/unsimplified step in the proof, no error in the Lemma. Confirmed numerically below.
2. d(alpha) (E-196, second item, l. 1124 and 1140). No extra hypothesis is needed. Elements of A are arrows of degree 0 (l. 321: "A the collection of all arrows of degree 0 which start in i"), the grading is non-negative (Thm l. 262), and d lowers degree by 1, so d(alpha) lies in degree -1, i.e. d(alpha) = 0. Hence d((dφ/α)α) = d(dφ/α) α with no second term, which is what l. 1140 uses, and d(red dφ) = -Σ d(dφ/α) α lies in the span of paths starting in A, which is what l. 1124 asserts. Without dα = 0 both lines would fail (hand example: dφ = γα + ζε with dα = δε, dζ proportional to γδ gives ∂²φ != 0); it is the non-negative grading that rescues them.
3. What Sect. 8 needs. Used and visible in the text: (a) non-negative grading, (b) arrows of A are degree 0 and not loops (l. 321) (so t(α) != i and a path q = dφ/α does not start in A, needed for Δ(dφ/α) to make sense), (c) free algebra, d determined by Leibniz on arrows, d² = 0 on Q. Nothing about cycles, J, finiteness of dim, or relations. So "no degree-0 loop at i" plus the standing setup (l. 278-296) is all that Sect. 8 (the ∂² = 0 section) uses; E-189's wording stands, and the statement should also list "non-negatively graded" since dα = 0 comes from it.
4. Definition cross-check. The last line of Def. 3.2 (l. 390), "- alpha dec(dφ/β)" with no (-1)^c_i factor, is consistent: c_i(t(αφ)) = 0 gives -(−1)^0 = −. Replacing it by "+" breaks ∂²=0 (control).

Would refute: a random instance of the test class where ∂² != 0 with the printed signs; none in 7 082 letters.

## Evidence

`scholar_oppermann_sign.py` builds a random dg quiver (kQ with d on degree-1 arrows chosen at random in degree-0 paths, degree-2 arrows with dω = dℓ1·ℓ2 − ℓ1·dℓ2; d²_kQ = 0 checked on every arrow) and the letters of M (φ, αφ, φβ^-1, αφβ^-1, α^*), with ∂ on letters exactly as Def. 3.2 and Leibniz with ||.||-signs, dec/Δ expanded at every junction at i.

| last term of ∂(αφβ^-1) | d²_kQ fails | ∂² != 0 on letters | letters |
|---|---|---|---|
| printed: - α dec(dφ/β) | 0 | 0 | 7 082 |
| control: + α dec(dφ/β) | 0 | 788 | 7 082 |

`scholar_oppermann_lemma.py` tests the Lemma's second formula on random paths p (length 1-3, start at i, not starting in A) and each α in A: lhs ∂(Δp α^-1) computed by Leibniz on the expanded M-element, rhs = the statement.

| coefficient s in "s (-1)^{c_i(t p)} Δ(dp/α)" | agree | disagree |
|---|---|---|
| s = -1 (statement l. 1058) | 610 | 0 |
| s = +1 (control) | 330 | 280 |

(610 = 337 end-at-i with start loop + 16 end-at-i + 182 + 75 end elsewhere; includes the multi-letter case with p_1 a loop at i, the one l. 1096-1105 treats.) The test also reproduces ∂² = 0 for φα^-1 and for degree-1 and degree-2 arrows φ.

## Reproduction

```
.venv/bin/python workshop/rounds/069/scholar_oppermann_sign.py 200   # <1 s
.venv/bin/python workshop/rounds/069/scholar_oppermann_lemma.py      # <1 s
```

## Prior record

E-196 flags exactly the two observations; E-189 gives the hypothesis list. New: both are resolved (sign: proof typo, Lemma fine; dα = 0 by degree), and an executable d² test of Def. 3.2 exists. Grep of RETRACTIONS.md for 1504.02617: nothing. Line numbers: Lemma l. 1055, sign l. 1058, proof end l. 1105, Prop. l. 1109, d(alpha) uses l. 1124, 1140; Def. 3.2 l. 383-391 (all of the v1 .tex). The paper's own proof of ∂²=0 says the "αφβ^-1 case is the same calculation" (l. 1143); the test supports that for degrees <= 2.

## Code changed

None in the repo outside the two new scripts above (no tests touched).

## Next

- Citation text for E-189/E-196 may say: "Sect. 8 needs only: non-negative grading, A = degree-0 arrows out of i, none a loop; d(alpha) = 0 by degree" (checked against l. 262, 321).
- Remaining for T5 from the Oppermann side: Prop. 7.6 (cancellation, a sketch) and the quasi-isomorphism claim are unchanged and not covered here. theorist/skeptic: a cycle meeting D(i) is still the open case in the repo, not in the paper.
- Limit: my test covers degrees <= 2 and random d; it is evidence for the signs, not a proof of Thm 8.1.
