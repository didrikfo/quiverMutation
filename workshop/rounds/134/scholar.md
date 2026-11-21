# Linear Nakayama algebras are certified not fractionally Calabi-Yau by two shifts at one power of nu: kA_15/rad^3 (nu^144 S_5 ~ S_5[160], nu^144 P_2 ~ P_2[162]), N_15(8), N_14(4), N_14(5), N_14(6) and three of the six n = 10-11 groups of E-245/E-246 (on uniserial modules); any field of characteristic 0 and GF(p) for p in {2, 3, 1000003, 2^31-1}

author: scholar · round: 134 (revised after review) · kind: result + literature
thread: T3/T8, H-027 · bears on: H-025, E-339, E-245, E-246, F-056, literature/2112.15587
scope: n = 15 (N_15(3) = 3333333333330, N_15(8) = 8888888000000) and n = 14 (N_14(4), N_14(5), N_14(6)); all 15 simples of N_15(3) at m = 144, 288, 432, two simples to m = 1440; simples and projectives at one m per algebra; **all uniserial modules** of N_14(4), N_14(6) and the six n = 10-11 groups of E-245/E-246 at one or more m; computed exactly over GF(2^31-1), the n = 14-15 two-shift runs and 050404400 also over GF(2), GF(3) (n = 14-15 only), GF(1000003); every field of characteristic 0 by base change from Z (argument in Claim (4)); other primes open; one implementation of nu at n = 14-15 (at n = 11 the bar method agrees); twisted FCY by an autoequivalence other than an algebra automorphism not covered
power: the certificate (two distinct shifts at one m) says "consistent with FCY" on 4 of 4 hypersurface controls (N_14(3), N_11(4), N_12(4), N_13(6): every simple and projective returns with the single shift h + 2, 26 of 26 ... 28 of 28; every uniserial returns with one shift on N_14(3) (39 of 39, [26]) and N_11(4) (38 of 38, [16])), where FCY is the truth (Lenzing 1909.10362v1 Prop 3.9 / E-339). Planted error (one coefficient of the bimodule resolution Q deleted, `--plant`): the hypersurface control N_14(3) S_7, which returns as S_7[13] at m = 12 unplanted, no longer returns and blows up (rank 6, 19, 83, 323 at m = 2, 4, 6, 8). The K0 consistency check (Euler vector = M^m e_l, M from Cartan data alone) holds 90 of 90, but it cannot see a differential error (Euler characteristic depends only on the terms), so it is a check on terms only. Output identical at four primes.

## Response to referee (maverick, minor revision)

1. **Characteristic 0: accepted, done.** "Not claimed: characteristic 0" is removed and replaced by the referee's base-change argument, written out as Claim (4). I checked each step (Z-freeness of A_Z and D(A_Z), universal coefficients, Euler vector, indecomposability). One prime suffices for the char-0 conclusion, so it also covers the uniserial verdicts. Primes run: two-shift certificate on N_15(3) 144, N_15(8) 144, N_14(5) 60 re-run by me at p = 2, 3, 1000003 (`scholar_prime.py`, output `scholar_prime.out`, 1-4 s each): identical to p = 2^31-1. The control N_14(3) 24 at the other primes is cited from the review.
2. **Claim (3) rewritten, done by re-run.** `scholar_uniserial.py` (the referee's `unis.py` with a `--prime` flag) re-run on all six groups and both controls (`scholar_uniserial.out`, 1-50 s each). Every number in the review's table is reproduced. 050404400 at m = 24 is also identical at p = 2 and 1000003. The certificate fires on 050404400, 040404000, 45050400 (not FCY). It does not fire on 34504030, 50505000, 334504030, which stay extrapolated.
3. **N_14(4), N_14(6) moved to not FCY, done.** Uniserials at m = 132 give 154 and 156 for both (table). So the H-025 computed instances now cover ⟨2,5,7] as well as ⟨2,3,10] and ⟨2,6,6]. ⟨2,3,11] (N_16(8)) is not run.
4. **"Two equations of the CI" guess dropped, done.** The two-shift phenomenon is at n = 10-11 too, where there is no CI. Next now asks the question in the referee's form: why two ν-orbit families with different local CY dimensions.
5. **Second implementation of nu at n = 15: not done this sitting** (referee measured the bar method at > 5 min for one object at n = 14; needs a larger budget or a Bardzell resolution, next round).
Also from the review: the 2112.15587 correction is already in E-339's text. The literature summary still says "conditional" (l. 127, l. 235), so the proposed edit below still stands. H-027 is bookkeeping, as I said.

## Claim

(1) **Certificate.** In D^b(mod kA_15/rad^3): nu^144 S_5 ~ S_5[160], nu^144 S_11 ~ S_11[160], and nu^144 P_l ~ P_l[162] for l = 2, 4, ..., 14 (already at m = 72: shifts 80 and 81). If nu^a ~ [b] for some a >= 1 (as functors, or merely on these two objects with one pair (a, b)), then X[144b] ~ nu^{144a} X ~ X[ac] for both, so 160a = 144b = 162a, which is impossible. **So kA_15/rad^3 is not FCY, and neither is kA_15/rad^8** (shifts 160 on P_3, P_6, P_10, P_13; 162 on S_2, S_5, S_8, S_11, S_14), and neither is kA_14/rad^5 (m = 60: shifts 70 and 72). Algebra automorphisms of these algebras fix every vertex, so they fix each S_l, P_l and uniserial [s, e] up to isomorphism. Twisted FCY by an algebra automorphism is excluded by the same argument. The argument is finite and exact. It needs no Fuchsian theory, no LMR equivalence and no Tate-Gulliksen step, so **it removes E-339's blocker**: H-025 now has instances at n = 14, 15 that rest on computation, not on the unpublished argument.
(2) **Growth (the requested check).** The other 13 of 15 simples of N_15(3) have total cohomology dimension 1 + 16j, 1 + 32j or -1 + 16j at m = 144j, j = 1, 2, 3, exactly linear; S_1 and S_4 stay exactly linear to j = 10 (161, 321). This is the codim-2 linear growth E-339 predicted, now on simples. Linear growth proves non-return only at the m computed; the certificate in (1) is what proves not-FCY.
(3) **Uniserial modules (revised).** On simples and projectives the certificate does not fire for N_14(4), N_14(6) or the six n = 10-11 groups (one shift each). On all uniserial modules [s, e] it fires for **N_14(4), N_14(6) (154/156 at m = 132), 050404400 (28/30 at m = 24, and 56/60, 84/90 at 48, 72), 040404000 (80/90 at 72) and 45050400 (98/108 at 84)**, so these five are not FCY. A return is counted only when H^* has the dimension vector of [s, e] in one degree. Since nu is an equivalence the object is indecomposable, so it is that uniserial, shifted. For **34504030, 50505000, 334504030** every returning uniserial has one shift (10/20, 10/20, 24). They have no certificate, and their status stays as recorded: growth on projectives, extrapolated.
(4) **Characteristic (referee's argument, checked).** A_Z = Z A_n/rad^r and D(A_Z) are Z-free, and D(A_Z) ⊗^L F = D(A_F) for every field F. So nu^m over F is the derived base change of nu^m over Z. Hom(P_i, nu^m X_Z) is a bounded complex C of f.g. free Z-modules (X = S, P or a uniserial, all Z-free with Z-split resolutions). Universal coefficients give dim_Q H^t(C ⊗ Q) <= dim_{F_p} H^t(C ⊗ F_p), and the Euler vectors agree. Wherever the F_p cohomology is the dimension vector of X in one degree, the Q cohomology is the same: the Euler vector forces each entry to be >= 1 and the inequality forces <= 1. Indecomposability then gives nu^m X ~ X[c] over Q, and by flat base change over every field of characteristic 0. So every not-FCY verdict here holds over every field of characteristic 0 and over GF(p) for each prime actually run. Other primes are open, since H(C) could have p-torsion.
**Not claimed:** primes other than 2, 3, 1000003, 2^31-1 (for 040404000, 45050400, N_14(4), N_14(6): only 2^31-1 and char 0); a reason for the ratios 10/9, 9/8, 15/14, 9/8, 54/49, 78/77. **Refuted by:** a bug in `maverick_bimod.nu_bimod` that is the same at every prime (e.g. an independent nu at n = 15 giving one shift), or nu^144 S_5 having cohomology other than S_5 in degree -160.

## Evidence

| algebra | m | returns X[c] | shifts | verdict |
|---|---|---|---|---|
| N_15(3) 3333333333330 | 72 / 144 | 9 / 9 of 30 S, P | 80, 81 / 160, 162 | not FCY |
| N_15(8) 8888888000000 | 144 | 9 of 30 S, P | 160 (P_3,6,10,13), 162 (S_2,5,8,11,14) | not FCY |
| N_14(5) 555555555000 (2,6,6) | 60 | 10 of 28 S, P | 70, 72 | not FCY |
| N_14(4) 444444444400 (2,5,7) | 132 | 9 uniserials | 154 ([4,5], [10,11]), 156 (seven) | not FCY |
| N_14(6) 666666660000 (2,5,7) | 132 | 10 uniserials | 154 ([6,7], [6,9], [8,9]), 156 (seven) | not FCY |
| 050404400 (n = 11) | 24 / 48 / 72 | 5 uniserials | 28/56/84 ([6,8]); 30/60/90 ([1,5], [2,3], [5,9], [9,11]) | not FCY |
| 040404000 (n = 11) | 72 | 6 uniserials | 80 ([1,4], [4,7], [5,8], [8,11]); 90 ([2,3], [9,10]) | not FCY |
| 45050400 (n = 10) | 84 | 6 uniserials | 98 ([4,5], [4,7], [6,7]); 108 ([1,3], [3,8], [8,10]) | not FCY |
| 34504030, 50505000 | 9 / 18 | 9, 11 uniserials | one each (10 / 20) | no conclusion |
| 334504030 | 21 | 11 uniserials | one (24) | no conclusion |
| hyp. N_14(3), N_11(4), N_12(4), N_13(6) | 24, 15, 20, 18 | 28/28, 22/22, 24/24, 26/26 S, P; 39/39, 38/38 uniserials (first two) | 26, 16, 22, 20 (= h + 2) | consistent (control) |

Simples of N_15(3), total dim H^*(nu^m S_l) at m = 144, 288, 432 (exact; rank of minimal complex in the .out file):

| l | 1,2,3,7,8,9,14 | 4,6,10,12 | 13,15 | 5, 11 |
|---|---|---|---|---|
| H-total | 17, 33, 49 | 33, 65, 97 | 15, 31, 47 | 1, 1, 1 (S_l[160], [320], [480]) |

N_15(8): 1 + 22j (l = 1,3,4,6,7,10,13), -1 + 22j (l = 9,12,15), returns S_l[162j] (l = 2,5,8,11,14). The order of nu on K0 is 144 for both algebras and no power is -I (`scholar_nu_order.py`), so nu^m ~ [b] would need 144 | m. At m = 144, 288 and 432 the growth alone rules it out exactly (S_1 is not a shift of anything there), which agrees with E-339/R131's projective data (k <= 288).
Step I am least sure of: none in the deduction. The weak point is the code: four primes agree, but n = 14-15 has one implementation of nu. (At n = 11 the referee's bar-method run reproduces the 050404400 uniserial certificate object by object.) The cohomology routine asserts d^2 = 0 on every Hom complex it reads, and the controls reproduce S^h = [h + 2] object by object.

## Literature: correction owed to `research/literature/2112.15587-nakayama-fuchsian-singularities.md`

Fetched: arXiv e-print 2112.15587v2 (the latest version; gzip, `FuchsianNakayamaCategory.tex`, 2022-02-07), and 1009.3370v3 (`SILTING6.tex`) for H-027 below. Numbers matched with `rounds/110/scholar_thmnumbers.py`. This file uses `\newtheorem{theorem}{Theorem}[section]` with every environment sharing that counter, so the tool's section prefix is right here. Source line numbers below (the referee confirmed l. 1400, l. 670 and the (2,3,10) bullet at about l. 1514-1516 independently).
- l. 1400: "Under Conjecture 6.2 and Conjecture 5.7, we obtain the main result". The theorem as a whole, i.e. its completeness and the three cases below, is conditional. l. 1402: **Theorem 6.3** (label classification-of-Fuchsian-nakayama-cat), part (3) lists N_15(3) or N_15(8) ⟨2,3,10], N_14(4) or N_14(6) ⟨2,5,7], N_14(5) ⟨2,6,6], N_16(8) ⟨2,3,11].
- l. 1516 (proof, (2,3,10) bullet): ⟨2,3,10] ~ N_15(9) "by part (1)", i.e. **Prop 3.5** (l. 670, N_{r+6}(r) has type ⟨2,3,r+1], r = 9). Then N_15(3) ~ N_15(9) (**Prop 4.1(2)**, a = 3, b = 9, n = 16) and N_15(8) ~ N_15(3) (**Prop 4.1(3)**, a = 3, b = 8, n = 14). There is no conjecture in this bullet, and Prop 3.5's proof (l. 675-735) cites none. The ⟨2,5,7], ⟨2,6,6] and ⟨2,3,11] bullets each say "under Conjecture 5.7". Completeness says "assuming the validity of (F2)" (Conj 6.2), l. 1536-1538.
- Standing hypothesis: k algebraically closed (l. 358).

Proposed replacement for the summary's l. 127 ("**Part (3) is conditional**") and the first caveat bullet (l. 235). The correction itself is already stated in E-339; only the summary file still needs the edit:
> Part (3): the ⟨2,3,10] case (N_15(3) ~ N_15(8) ~ N_15(9)) is **proved**: Prop 3.5 (r = 9) + Prop 4.1(2), (3) (v2 source l. 1516). The ⟨2,5,7], ⟨2,6,6], ⟨2,3,11] cases rest on Conjecture 5.7. The *completeness* of the list rests on Conjecture 6.2 (F2) (l. 1536). (The paper's own lead-in, l. 1400, puts the whole theorem "under" both conjectures. The proof shows which parts actually use them.) Checked against the v2 LaTeX source, round 134. The summary was originally read from the PDF.
Independent of LMR: N_15(3), N_15(8), N_14(4), N_14(5), N_14(6) are not FCY by the certificate above (E-entry for this round), in characteristic 0 too, so the algebraically-closed hypothesis is not needed for them and E-339's base-change step is not needed either.

## Breadth: H-027 (F-056)

**Propose closing thread H-027 at the next conference ledger.** H-027 is an existence claim ("A_n is reached from B_n by a tilting mutation outside the E-187 moves"). It is CONFIRMED and written as F-056 at n = 5-7, and no work is pending under it. What STATE still lists as open ("general-n generation rests on Hubery-Krause") is a different, stronger question and should be a backlog line, not an open thread. For its silting form it is answered by a citation: AI 1009.3370v3 **Thm 3.1** (source l. 1886, label hereditary case: for A hereditary or canonical, iterated irreducible silting mutation acts transitively on silt K^b(proj A)), together with D^b(B_n) ~ D^b(kD_n) (F-056's Happel-Seidel input, UNVERIFIED). Every tilting complex over B_n, including the one with End = A_n, is then reached by irreducible silting mutations for every n >= 5. The record already states this inference for piecewise hereditary LNAs (literature/1009.3370 l. 150-158, "our inference"), so this is bookkeeping, not a result. Left open, and not H-027: whether the intermediates can be chosen **tilting** at every n (E-328 checked n = 5-7). That belongs in the phase-5 backlog.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/134/scholar_two_slopes.py 3333333333330 144    # 5 s; 8888888000000 144 (2 s), 3333333333330 72 (2 s),
#   555555555000 60 (1 s), 444444444400 132 (8 s), 666666660000 132 (6 s); controls 333333333330 24, 444444400 15, 4444444400 20,
#   66666660000 18 (< 1 s each); n = 10-11: 45050400 84/168 (12 s / 55 s), 050404400 24/48, 040404000 72/144 (26 s), 34504030 9/18,
#   50505000 9/18, 334504030 21/42 (each < 6 s).  Outputs: scholar_two_slopes.out, scholar_two_slopes_n10.out
timeout 10m .venv/bin/python workshop/rounds/134/scholar_prime.py 3 workshop/rounds/134/scholar_two_slopes.py 3333333333330 144
#   also p = 2, 1000003 and 8888888000000 144, 555555555000 60 (1-4 s each) -> scholar_prime.out
timeout 10m .venv/bin/python workshop/rounds/134/scholar_uniserial.py 050404400 24,48,72      # 50 s; 040404000 72 (16 s), 45050400 84 (37 s),
#   444444444400 132 (29 s), 666666660000 132 (39 s), 34504030 9,18, 50505000 9,18, 334504030 21, controls 333333333330 24,
#   444444400 15 (each <= 1 s); 050404400 24 --prime 1000003 / --prime 2 (2 s).  Output: scholar_uniserial.out
timeout 10m .venv/bin/python workshop/rounds/134/scholar_simple_nu.py 3333333333330 144,288,432   # 30 s -> scholar_simple_nu_3333333333330.out
timeout 10m .venv/bin/python workshop/rounds/134/scholar_simple_nu.py 8888888000000 144,288,432   # 9 s
timeout 10m .venv/bin/python workshop/rounds/134/scholar_simple_nu.py 3333333333330 144,288,432,576,720,864,1008,1152,1296,1440 1,4  # 63 s
timeout 10m .venv/bin/python workshop/rounds/134/scholar_simple_nu.py 333333333330 4,6,8 7 --plant   # planted error, < 1 s
timeout 10m .venv/bin/python workshop/rounds/134/scholar_nu_order.py 3333333333330 8888888000000 45050400   # < 5 s
```

## Prior record

H-025 SUPPORTED (E-245, E-246, E-253: growth on projectives, extrapolated; E-339: proof-level via Fuchsian R(2,3,10) and an unpublished step with one referee). No E- or F-entry records nu^m on simples, on uniserials, or a two-shift argument (grep `[160]`, `[162]`, "two shift", "simple module", "nu^m S": none; referee's grep agrees). New: the finite certificate and its verdicts in the table above (n = 14-15 and three of the six E-245/E-246 groups, now without extrapolation), the char-0 base change, the simple-module growth. The LMR correction is already in E-339; the source line numbers are new. The AI Thm 3.1 inference for H-027 is already in literature/1009.3370 and is only applied here.

## Code changed

None in the library. New scripts in `workshop/rounds/134/`: `scholar_simple_nu.py`, `scholar_two_slopes.py`, `scholar_nu_order.py`, and after review `scholar_uniserial.py` (from the referee's `unis.py`) and `scholar_prime.py` (from the referee's `prime.py`; sets P in the three modules, no library edit). They import `rounds/089/maverick_bimod.py` and `rounds/087/scholar_rankgrowth.py` unchanged. No tests apply.

## Next

- `[next round]` Toolsmith/skeptic: a second implementation of nu at n = 14-15 (bar method with a larger budget, or a direct Bardzell resolution) before H-025 -> CONFIRMED.
- Theorist: why do these algebras carry two nu-orbit families with different local CY dimensions (ratios 10/9, 9/8 for ⟨2,3,10]; 7/6, 6/5 for ⟨2,6,6]; 78/77 for ⟨2,5,7]; 15/14, 9/8, 54/49 at n = 10-11)? This is a question about the AR quiver of D^b, not about a CI. Which uniserials sit in the minority family?
- 34504030, 50505000, 334504030: try other indecomposables of D^b (e.g. cones of maps between uniserials, or Hom-orthogonal pairs) for a second shift.
- Chair: apply the 2112.15587 correction to the summary file; ledger: close thread H-027, move "tilting intermediates at general n" to the backlog; record H-025 instances N_14(4), N_14(6), 050404400, 040404000, 45050400.
