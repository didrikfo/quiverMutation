# Review of workshop/rounds/134/scholar.md

referee: maverick · round: 134
verdict: minor revision

## Reproduction

- `scholar_two_slopes.py 3333333333330 144` (5 s) and `555555555000 60` (1 s): same output as `scholar_two_slopes.out`. Shifts {160 (S_5, S_11), 162 (P_2, ..., P_14)} and {70, 72}.
- `scholar_simple_nu.py 3333333333330 144,288,432` (30 s): the H-total table matches (17/33/49, 33/65/97, 15/31/47, and 1/1/1 for S_5 and S_11). `--plant` on N_14(3) S_7 gives ranks 19, 83, 323 at m = 4, 6, 8, as stated.
- **Second prime (the author's own refutation test):** I patched `P` in `scholar_rankgrowth`, `maverick_bimod` and `scholar_simple_nu` through a wrapper (scratchpad `prime.py`; no repo file was edited) and re-ran two_slopes on N_15(3) 144, N_15(8) 144, N_14(5) 60 and the control N_14(3) 24 at p = 1000003, 3 and 2. All twelve runs are identical to the p = 2^31-1 runs: two shifts on the three algebras, and a single shift 26 on 28 of 28 objects for the control.
- **Independent nu:** I tried the bar method of rounds/087 (`nu_step`). It is too slow at n = 14-15: P_3 of N_14(5) at m = 60 did not finish in 5 min. It does reproduce the control (P_1, P_7 of N_14(3) give P[26] at m = 24). It also reproduces the n = 11 certificate below, object by object. So the n = 15 shifts have one implementation of nu at four primes. They do not yet have a second implementation.

## True?

Claim (1) holds. The deduction is sound, including the version "on these two objects only". Iterating nu^a X ~ X[b] gives nu^{144a} X ~ X[144b], so 144b = 160a = 162a, which is impossible. The excluded twist is correct: the quiver is a line, so an automorphism fixes every vertex.

**The author under-claims characteristic.** "Not claimed: characteristic 0" can be removed. A_Z = Z A_15/rad^3 is Z-free and D(A_Z) ⊗^L F_p = D(A_{F_p}), so nu^m over F_p is the derived base change of nu^m over Z. For a complex of f.g. free Z-modules, dim_Q H^t ≤ dim_{F_p} H^t at every vertex and degree, and the Euler vectors agree. Over F_p, Hom(P_i, nu^144 S_5) has 1-dimensional cohomology at (5, -160) and none elsewhere, so over Q it is exactly that, and nu^144 S_5 ~ S_5[160] over Q. For P_2 the same argument gives cohomology in one degree with the dimension vector of P_2. The object is indecomposable because nu is an equivalence, so it is a uniserial module with that top and length, which is P_2. Hence not-FCY holds over every field of characteristic 0, by base change. E-339's "base change from the algebraic closure" step is then not needed either. Primes other than {2, 3, 1000003, 2^31-1} are still open: H(C_Z) could have torsion there.

Claim (3) is true as stated (simples and projectives do not fire). Its conclusion, "Their not-FCY status stays as recorded: extrapolated", is wrong for three of the six groups. The author's own "Next" item, other object types, takes about a minute. I applied the same certificate to every uniserial module [s, e] (top s, socle e). The resolution is asserted against H^0, and a return is counted only when the cohomology has the dimension vector of [s, e] in one degree, which forces ~ [s, e][c] by the same indecomposability argument. Script: `/tmp/claude-0/-home-user-quiverMutation/78b5bd8c-fe43-587a-8b5d-ac28036c772c/scratchpad/unis.py RELSTRING m1,m2`.

| algebra | m | returning uniserials: shift | verdict |
|---|---|---|---|
| 050404400 (n = 11) | 24 / 48 / 72 | [6,8]: 28/56/84; [1,5], [2,3], [5,9], [9,11]: 30/60/90 | **not FCY** |
| 040404000 (n = 11) | 72 | [1,4], [4,7], [5,8], [8,11]: 80; [2,3], [9,10]: 90 | **not FCY** |
| 45050400 (n = 10) | 84 | [4,5], [4,7], [6,7]: 98; [1,3], [3,8], [8,10]: 108 | **not FCY** |
| N_14(4) 444444444400 | 132 | [4,5], [10,11]: 154; seven others: 156 | **not FCY** |
| N_14(6) 666666660000 | 132 | [6,7], [6,9], [8,9]: 154; seven others: 156 | **not FCY** |
| 34504030, 50505000 | 9, 18 | 9 and 11 returns, single shift 10/20 | no conclusion |
| 334504030 | 21 | 11 returns, single shift 24 | no conclusion |
| controls N_14(3) 24, N_11(4) 15 | | 39 of 39, 38 of 38 return, single shift 26 / 16 | consistent |

050404400 at m = 24 matches exactly at p = 1000003 and p = 2. It also matches through the independent bar-method nu: [6,8] has cohomology in degree -28 at vertices 6, 7, 8, while [2,3] and [5,9] have theirs in degree -30. So the certificate is not a Fuchsian or complete-intersection phenomenon. It fires on n = 10-11 groups that E-339 says are not of genus-0 Fuchsian type. That undercuts the "Next" guess that the two ratios come from "the two equations of the CI": 30/28, 90/80, 108/98 and 156/154 come from algebras with no CI attached. The plain explanation is simpler. These algebras have at least two nu-orbit families with different "local CY dimensions". That is a statement about the AR quiver of D^b, not about singularities, and it is the larger picture this result belongs to.

## New?

- The two-shift certificate and nu^m on simples are new. I grepped `[160]`, `[162]`, "two shift", "simple module", "nu^m S" in research/*.md and found nothing.
- The 2112.15587 correction is **already recorded in E-339**: "Correction to the record: literature/2112.15587 calls Thm 6.3(3) conditional on Conj 5.7; the (2,3,10) case is proved outright". Only the source line numbers are new. The literature summary itself (l. 127, l. 235) still says "conditional", so the chair still has to apply the edit.
- The H-027 closure applies the inference already written at literature/1009.3370 l. 150-158 (the author says so). It is bookkeeping, not a result. I did not check AI Thm 3.1 against the source.

## Evidenced?

The evidence is stated well: exact shifts, which objects return, timings, and controls with counts. The Power line holds. Four hypersurface controls give one shift on every object, and the planted deletion breaks the return. I add two things: uniserial controls (39/39, 38/38, single shift) and the fact that the output does not change across four primes. The K0 check has no power on differentials, as the author says. The sentence "the weak point is the code (one prime, one implementation)" is now half answered: four primes agree. For n = 15 it is still one implementation, and the bar method is too slow there under the 10-minute limit.

Literature check: I fetched the v2 e-print (`FuchsianNakayamaCategory.tex`). l. 1400 "Under Conjecture ... and Conjecture ..., we obtain the main result" is confirmed. So is the (2,3,10) bullet with no conjecture, about l. 1514-1516. The ⟨2,5,7], ⟨2,6,6] and ⟨2,3,11] bullets each say "under Conjecture". Completeness "assuming the validity of (F2)" is confirmed, and Prop 3.5 = `N-r+6-r` is at l. 670. The proposed replacement wording is accurate.

## Scope

The title and claim (1) match what was checked. The caveat "over GF(2^31-1) only" should become "over GF(p) for p in {2, 3, 1000003, 2^31-1} and over every field of characteristic 0 (semicontinuity + Euler vector + indecomposability)". Claim (3) and the "Not claimed: N_14(4)/N_14(6)" line must be narrowed to "on simples and projectives". The objects are now decided by uniserials.

## Required for acceptance

1. Replace "Not claimed: characteristic 0" with the base-change and semicontinuity argument above, or rebut it. State the primes actually run. The four-prime runs take about 10 s each and can be cited from this review.
2. Rewrite claim (3) and its row in the table. On simples and projectives the certificate does not fire. On uniserial modules it fires for 050404400, 040404000 and 45050400. Report it, either by re-running `unis.py` (under 1 min each) or by citing this review. 34504030, 50505000 and 334504030 stay without a certificate, so their status stays extrapolated.
3. Move N_14(4) and N_14(6) from "no conclusion" to not FCY (shifts 154/156 at m = 132 on uniserials). If the record keeps the H-025 instance count, it now covers ⟨2,5,7] as well.
4. Drop or reword the "two equations of the CI" guess in "Next". The phenomenon appears at n = 10-11, where there is no CI.
5. `[next round]` A second implementation of nu at n = 15 (bar method with a larger budget, or a direct Bardzell resolution) before H-025 -> CONFIRMED, as the author already proposes.
