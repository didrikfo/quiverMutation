# kA_15/rad^3 is not fractionally Calabi-Yau, by a finite exact certificate at n = 15 (and n = 14, rad^5): nu^144 S_5 ~ S_5[160] while nu^144 P_2 ~ P_2[162]; 13 of 15 simples grow linearly in m

author: scholar · round: 134 · kind: result + literature
thread: T3/T8, H-027 · bears on: H-025, E-339, E-245, E-246, F-056, literature/2112.15587
scope: n = 15 (N_15(3) = 3333333333330, N_15(8) = 8888888000000) and n = 14 (N_14(5)); all 15 simples at m = 144, 288, 432, two simples to m = 1440; simples and projectives at one m per algebra; exact over GF(2^31-1) only (no second prime, no char 0); E-245/E-246's six n = 10-11 groups tested with the same certificate (it does not fire there); twisted FCY by an autoequivalence other than an algebra automorphism not covered
power: the certificate (two distinct shifts at one m) says "consistent with FCY" on 4 of 4 hypersurface controls (N_14(3), N_11(4), N_12(4), N_13(6): every simple and projective returns with the single shift h + 2, 26 of 26 ... 28 of 28), where FCY is the truth (Lenzing 1909.10362v1 Prop 3.9 / E-339). Planted error (one coefficient of the bimodule resolution Q deleted, `--plant`): the hypersurface control N_14(3) S_7, which returns as S_7[13] at m = 12 unplanted, no longer returns and blows up (rank 6, 19, 83, 323 at m = 2, 4, 6, 8). The K0 consistency check (Euler vector = M^m e_l, M from Cartan data alone) holds 90 of 90, but it cannot see a differential error (Euler characteristic depends only on the terms), so it is a check on terms only.

## Claim

(1) **Certificate.** Over k = GF(2^31-1), in D^b(mod kA_15/rad^3): nu^144 S_5 ~ S_5[160], nu^144 S_11 ~ S_11[160], and nu^144 P_l ~ P_l[162] for l = 2, 4, ..., 14 (already at m = 72: shifts 80 and 81). If nu^a ~ [b] for some a >= 1 (as functors, or merely on these two objects with one pair (a, b)), then X[144b] ~ nu^{144a} X ~ X[ac] for both, so 160a = 144b = 162a, which is impossible. **So kA_15/rad^3 is not FCY over this field, and neither is kA_15/rad^8** (shifts 160 on P_3, P_6, P_10, P_13; 162 on S_2, S_5, S_8, S_11, S_14), and neither is kA_14/rad^5 (m = 60: shifts 70 and 72). Algebra automorphisms of these algebras fix every vertex, so they fix each S_l and P_l up to isomorphism and twisted FCY by an algebra automorphism is excluded by the same argument. The argument is finite and exact. It needs no Fuchsian theory, no LMR equivalence and no Tate-Gulliksen step, so **it removes E-339's blocker**: H-025 now has an instance at n = 15 that rests on computation (two derived-equivalent algebras, both certified), not on the unpublished argument.
(2) **Growth (the requested check).** The other 13 of 15 simples of N_15(3) have total cohomology dimension 1 + 16j, 1 + 32j or -1 + 16j at m = 144j, j = 1, 2, 3, exactly linear; S_1 and S_4 stay exactly linear to j = 10 (161, 321). This is the codim-2 linear growth E-339 predicted, now on simples. Linear growth proves non-return only at the m computed; the certificate in (1) is what proves not-FCY.
(3) **Null.** On the six n = 10-11 groups of E-245/E-246 (45050400, 050404400, 040404000, 34504030, 50505000, 334504030) the returning simples and projectives have **one** shift at m = order and 2 x order, so the certificate does not fire. Their not-FCY status stays as recorded: growth on projectives, extrapolated.
**Not claimed:** characteristic 0 or other primes (the computation is mod p; the integer complexes could in principle have different ranks at another prime); anything about N_14(4)/N_14(6) (single shift 156 at m = 132); a reason for the ratios 10/9 and 9/8. **Refuted by:** a bug in `maverick_bimod.nu_bimod`, for example a run at a second prime giving one shift, or nu^144 S_5 having cohomology other than S_5 in degree -160.

## Evidence

| algebra | m | returns X[c] (of 2n) | shifts | verdict |
|---|---|---|---|---|
| N_15(3) 3333333333330 | 72 / 144 | 9 / 9 | 80, 81 / 160, 162 | not FCY |
| N_15(8) 8888888000000 | 144 | 9 | 160 (P_3,6,10,13), 162 (S_2,5,8,11,14) | not FCY |
| N_14(5) 555555555000 (2,6,6) | 60 | 10 | 70, 72 | not FCY |
| N_14(4) / N_14(6) | 132 | 6 / 6 | 156 / 156 | no conclusion |
| hyp. N_14(3), N_11(4), N_12(4), N_13(6) | 24, 15, 20, 18 | 28/28, 22/22, 24/24, 26/26 | 26, 16, 22, 20 (= h + 2) | consistent (control) |
| 45050400, 050404400, 040404000 | order, 2 x order (84, 24, 72) | 2 each | one each (108/216, 30/60, 80/160) | no conclusion |
| 34504030, 50505000, 334504030 | 9/18, 9/18, 21/42 | 4, 4, 6 | one each (10/20, 10/20, 24/48) | no conclusion |

Simples of N_15(3), total dim H^*(nu^m S_l) at m = 144, 288, 432 (exact; rank of minimal complex in the .out file):

| l | 1,2,3,7,8,9,14 | 4,6,10,12 | 13,15 | 5, 11 |
|---|---|---|---|---|
| H-total | 17, 33, 49 | 33, 65, 97 | 15, 31, 47 | 1, 1, 1 (S_l[160], [320], [480]) |

N_15(8): 1 + 22j (l = 1,3,4,6,7,10,13), -1 + 22j (l = 9,12,15), returns S_l[162j] (l = 2,5,8,11,14). The order of nu on K0 is 144 for both algebras and no power is -I (`scholar_nu_order.py`), so nu^m ~ [b] would need 144 | m. At m = 144, 288 and 432 the growth alone rules it out exactly (S_1 is not a shift of anything there), which agrees with E-339/R131's projective data (k <= 288).
Step I am least sure of: none in the deduction. The weak point is the code (one prime, one implementation of nu). The cohomology routine asserts d^2 = 0 on every Hom complex it reads, and the hypersurface controls reproduce S^h = [h + 2] object by object.

## Literature: correction owed to `research/literature/2112.15587-nakayama-fuchsian-singularities.md`

Fetched: arXiv e-print 2112.15587v2 (the latest version; gzip, `FuchsianNakayamaCategory.tex`, 2022-02-07), and 1009.3370v3 (`SILTING6.tex`) for H-027 below. Numbers matched with `rounds/110/scholar_thmnumbers.py`. This file uses `\newtheorem{theorem}{Theorem}[section]` with every environment sharing that counter, so the tool's section prefix is right here. Source line numbers below.
- l. 1400: "Under Conjecture 6.2 and Conjecture 5.7, we obtain the main result". The theorem as a whole, i.e. its completeness and the three cases below, is conditional. l. 1402: **Theorem 6.3** (label classification-of-Fuchsian-nakayama-cat), part (3) lists N_15(3) or N_15(8) ⟨2,3,10], N_14(4) or N_14(6) ⟨2,5,7], N_14(5) ⟨2,6,6], N_16(8) ⟨2,3,11].
- l. 1516 (proof, (2,3,10) bullet): ⟨2,3,10] ~ N_15(9) "by part (1)", i.e. **Prop 3.5** (l. 670, N_{r+6}(r) has type ⟨2,3,r+1], r = 9). Then N_15(3) ~ N_15(9) (**Prop 4.1(2)**, a = 3, b = 9, n = 16) and N_15(8) ~ N_15(3) (**Prop 4.1(3)**, a = 3, b = 8, n = 14). There is no conjecture in this bullet, and Prop 3.5's proof (l. 675-735) cites none. The ⟨2,5,7], ⟨2,6,6] and ⟨2,3,11] bullets each say "under Conjecture 5.7". Completeness says "assuming the validity of (F2)" (Conj 6.2), l. 1536-1538.
- Standing hypothesis: k algebraically closed (l. 358).

Proposed replacement for the summary's l. 127 ("**Part (3) is conditional**") and the first caveat bullet (l. 235):
> Part (3): the ⟨2,3,10] case (N_15(3) ~ N_15(8) ~ N_15(9)) is **proved**: Prop 3.5 (r = 9) + Prop 4.1(2), (3) (v2 source l. 1516). The ⟨2,5,7], ⟨2,6,6], ⟨2,3,11] cases rest on Conjecture 5.7. The *completeness* of the list rests on Conjecture 6.2 (F2) (l. 1536). (The paper's own lead-in, l. 1400, puts the whole theorem "under" both conjectures. The proof shows which parts actually use them.) Checked against the v2 LaTeX source, round 134. The summary was originally read from the PDF.
Independent of LMR: N_15(3), N_15(8) and N_14(5) are not FCY by the certificate above (E-entry for this round).

## Breadth: H-027 (F-056)

**Propose closing thread H-027 at the next conference ledger.** H-027 is an existence claim ("A_n is reached from B_n by a tilting mutation outside the E-187 moves"). It is CONFIRMED and written as F-056 at n = 5-7, and no work is pending under it. What STATE still lists as open ("general-n generation rests on Hubery-Krause") is a different, stronger question and should be a backlog line, not an open thread. For its silting form it is answered by a citation: AI 1009.3370v3 **Thm 3.1** (source l. 1886, label hereditary case: for A hereditary or canonical, iterated irreducible silting mutation acts transitively on silt K^b(proj A)), together with D^b(B_n) ~ D^b(kD_n) (F-056's Happel-Seidel input, UNVERIFIED). Every tilting complex over B_n, including the one with End = A_n, is then reached by irreducible silting mutations for every n >= 5. The record already states this inference for piecewise hereditary LNAs (literature/1009.3370 l. 150-158, "our inference"). Left open, and not H-027: whether the intermediates can be chosen **tilting** at every n (E-328 checked n = 5-7). That belongs in the phase-5 backlog.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/134/scholar_two_slopes.py 3333333333330 144    # 5 s; 8888888000000 144 (2 s), 3333333333330 72 (2 s),
#   555555555000 60 (1 s), 444444444400 132 (8 s), 666666660000 132 (6 s); controls 333333333330 24, 444444400 15, 4444444400 20,
#   66666660000 18 (< 1 s each); n = 10-11: 45050400 84/168 (12 s / 55 s), 050404400 24/48, 040404000 72/144 (26 s), 34504030 9/18,
#   50505000 9/18, 334504030 21/42 (each < 6 s).  Outputs: scholar_two_slopes.out, scholar_two_slopes_n10.out
timeout 10m .venv/bin/python workshop/rounds/134/scholar_simple_nu.py 3333333333330 144,288,432   # 30 s -> scholar_simple_nu_3333333333330.out
timeout 10m .venv/bin/python workshop/rounds/134/scholar_simple_nu.py 8888888000000 144,288,432   # 9 s
timeout 10m .venv/bin/python workshop/rounds/134/scholar_simple_nu.py 3333333333330 144,288,432,576,720,864,1008,1152,1296,1440 1,4  # 63 s
timeout 10m .venv/bin/python workshop/rounds/134/scholar_simple_nu.py 333333333330 4,6,8 7 --plant   # planted error, < 1 s
timeout 10m .venv/bin/python workshop/rounds/134/scholar_nu_order.py 3333333333330 8888888000000 45050400   # < 5 s
```

## Prior record

H-025 SUPPORTED (E-245, E-246, E-253: growth on projectives, extrapolated; E-339: proof-level via Fuchsian R(2,3,10) and an unpublished step with one referee). No E- or F-entry records nu^m on simples or a two-shift argument (grep `[160]`, `[162]`, "simple module": none). New: the finite certificate and its verdicts in the table above, the simple-module growth, and the null on the n = 10-11 groups. The LMR correction was raised in R131 (E-339 scope line) and is now matched in source with line numbers. The AI Thm 3.1 inference for H-027 is already in literature/1009.3370 and is only applied here.

## Code changed

None in the library. New scripts in `workshop/rounds/134/`: `scholar_simple_nu.py`, `scholar_two_slopes.py`, `scholar_nu_order.py` (they import `rounds/089/maverick_bimod.py` and `rounds/087/scholar_rankgrowth.py` unchanged). No tests apply.

## Next

- Skeptic/toolsmith: rerun `scholar_two_slopes.py 3333333333330 144` at a second prime (P is a module constant in rounds/087), and through an independent nu (the bar method of rounds/087), before H-025 -> CONFIRMED at n = 15. Then a theorist referee can retire E-339's unpublished step as the load-bearing input.
- Theorist: why 10/9 and 9/8 for ⟨2,3,10] and 7/6, 6/5 for ⟨2,6,6]? (A guess, UNVERIFIED: two Gorenstein-parameter-like ratios from the two equations of the CI.) Why one shift for ⟨2,5,7]?
- n = 10-11 groups: look for two shifts among other objects (indecomposable modules, not only S and P); one more object type is cheap.
- Chair: apply the 2112.15587 correction; ledger: close thread H-027, move "tilting intermediates at general n" to the backlog.
