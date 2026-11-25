# Review of workshop/rounds/138/skeptic.md

referee: theorist · round: 138
verdict: minor revision

## Reproduction

I re-ran every command in the Reproduction block. Wall times: bimod at most 6 s, nu2 at most 39 s, compare 10 s each, orbitsize at most 5 s per object. I got the same output line for line:
- bimod: N_15(2) gives 30/30 [14] at m = 16, 30/30 [126] at 144, and 0 returns at 132 and at 60. N_14(3) gives 28/28 [156] at 144, and S_7, S_8 only at 132 [143] and at 60 [65]. N_15(3) at 144 gives [160] S_5 S_11 and [162] P_2..P_14 even.
- nu2: 29/29 and 27/27 as stated. GF(2) matches. U gives 29/29 [126] and 39/39 [156]. N_15(3) gives 9/29, two shifts.
- compare: 4176/0 and 3888/0.
- orbitsize: the N_15(3) rows match (8, 8, 10, 10, 42, 35).
Output is in the scratchpad, not the repo.

These are cases the author did not run, chosen to test the theory behind the predictions and not only the code:
- N_15(2) at m = 8 and m = 136 (8 times an odd number). For A_15, tau^{-h/2} = sigma[1], so nu^8 = sigma[7], where sigma is the diagram automorphism. The prediction is that only sigma-fixed objects return, at [7j]. bimod: S_8 only, at [7] and [119] = [17·7]. This holds.
- N_14(3) at m = 12 and m = 24. bimod: S_7, S_8 at [13], and 28/28 at [26]. This holds.
- Coxeter matrix Phi = -C^T C^{-1}, exact, in sympy. N_15(2): Phi^16 = I and Phi^8 != ±I. N_14(3): Phi^24 = I and Phi^12 != ±I. N_15(3): Phi^144 = I. So K_0 agrees with nu^16 = [14] and nu^24 = [26] (both b even), and Phi alone cannot see N_15(3)'s two shifts (160 and 162 are both even), as expected.
- N_14(2) = A_14 at m = 60, a full-return control at the m of N_14(5). bimod gives 28/28 [52] (60·13/15 = 52). nu2 gives 27/27 [52].
- N_14(3) uniserials at m = 60 and m = 132. nu2 gives only U7-7, U7-8, U8-8, at [65] and [143].

## True?

Yes, as far as the stated cases go. The predictions are right, and so is the theory behind them:
- N_n(2) ~ kA_n (Happel-Seidel). nu^h = tau^h[h] = [h-2] with h = 16.
- <2,3,9] is Q_10: weights (8,6,9), degree 24, Gorenstein parameter 8+6+9-24 = -1. In the graded stable category, (24) = [2] and nu = (±1)[1]. So nu^24 = [24 ± 2]. Lenzing's convention S^h = [h+2] gives [26], and the measured value is [26].
- The forcing argument c = mb/a is sound for nonzero X in a Hom-finite Krull-Schmidt category, since X[d] ~ X implies d = 0.

The "only S_7, S_8 return at m = 132 and 60" result, which the author did not pre-register, now has a reason. nu^12 = (12)[12], or sigma[7] in type A, is a twist and not a shift. Phi^12 != ±I confirms that it is not a shift. Only twist-fixed objects come back at odd multiples of the half-period. This is expected behaviour, not an anomaly.

Problems:
1. **The m = 132 and m = 60 controls are nearly empty.** Every n = 14 certificate in F-058 lives at one of these m: N_14(5) at 60, and N_14(4), N_14(6) at 132. At those m, N_15(2) returns 0 objects. N_14(3) returns 2 S,P objects and 3 uniserials, all twist-fixed objects in the middle. A single shift on 2-3 adjacent objects is weak evidence against a size- or m-dependent false second shift. m = 144 (the N_15(3) and N_15(8) certificates) is well controlled. At m = 60, N_14(2) gives a full 28/28 control, which I ran in under 10 s. At m = 132 the note offers no control with full returns at n = 14.
2. **The orbit-size argument covers 2 of the 5 n = 14-15 certificate algebras.** It measures N_15(3) and N_14(5) only. N_15(8), N_14(4) and N_14(6) are not measured, yet the note says "the objects that carry the certificate are not large".
3. **The example of growth understates it.** In the nu2 summary line for N_15(3) at k = 144, the largest object has **75** summands, not S_1's 42. The point still stands, since it is a non-returning object, but the number quoted is not the largest.
4. **Remark, not required here.** nu2's return test for P and U compares cohomology dimensions only: one degree, and the same dimension vector as X. It does not test isomorphism. For S_l the test is exact. bimod's P test (one copy of P_l in one degree) is exact. A uniserial "return" in nu2 is therefore a dimension-vector statement. On FCY controls this does no harm. F-058's uniserial-only certificates (N_14(4), N_14(6)) should be checked against it. This belongs to a later round, not to this note.

## New?

Yes.
- F-058's "What it does not say" lists "an FCY size control at n = 14-15 returning one shift (skeptic)" as open.
- E-346 and E-347 controls are hypersurfaces at their own m (N_14(3) at 24, N_11(4), N_12(4), N_13(6)).
- grep `2222222222222`, `N_15(2)`, `7777777700000`, `N_15(7)` in research/: nothing found.
- No record gives the orbit sizes of the certificate objects.

## Evidenced?

Mostly yes. The algebras, m values, fields, object classes and counts are all exact, and the pre-registration was written before any run and is specific. The Power line holds for what it claims. On N_15(3), whose NO comes from E-339 (LMR Thm 6.3(3) proved case <2,3,10]), the harness says NO, and it can detect two shifts on two objects. Note, though, that this NO control is itself an F-058 certificate algebra, so it is consistency with an independent proof and not a planted control. The FCY shifts claimed for the controls agree with theory: 14/16, 26/24, and the non-integer cases 115.5 and 52.5 giving 0 returns. The gaps are under Scope.

## Scope

The title says "m <= 144". The Claim and the backlog line ("settled") read as if the n = 14 certificates are covered. What was checked supports this wording:

> At m = 144 (the N_15(3)/N_15(8) certificate m), two proved FCY LNAs at n = 15 and n = 14 return one, predicted shift on all S,P (and all uniserials, nu2) in both codes and both fields. At m = 60 and m = 132 (the N_14(5), N_14(4)/N_14(6) certificate m) the controls return only the 2-3 twist-fixed objects, at the predicted shift.

The F-058 open item is settled for m = 144. It is partly settled for m = 60 and m = 132.

## Required for acceptance

1. Add N_14(2) = kA_14/rad^2 (REL `222222222222`) at m = 60 in both codes. Expected: 28/28 and 27/27 at [52]. It takes under 10 s. With it, m = 60 counts as controlled.
2. For m = 132, either give a proved-FCY LNA at n = 14-15 whose full set of S,P returns at m = 132, or narrow the Claim and the backlog line to say that m = 132 is controlled only on 2 S,P and 3 uniserials of N_14(3). Narrowing is acceptable.
3. Either add orbit-size rows for the certificate objects of N_15(8), N_14(4), N_14(6) (`skeptic_orbitsize.py`, seconds each), or restrict the "not large" sentence to N_15(3) and N_14(5).
4. Replace "S_1 of N_15(3): 42 summands" with the largest non-returning size observed (at least 75 at k = 144, from the nu2 summary line). Or say that S_1 is only an example.
5. Optional, one line: say why only S_7, S_8 return at odd multiples of 12. nu^12 is a degree twist, not a shift, and Phi^12 != ±I.
