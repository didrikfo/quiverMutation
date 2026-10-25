# nu^15 P_l = P_l[16] on projectives for the n = 11 group 535 (kA_11/rad^4, rad^5, 033450000); in the n = 10, 11 groups 34504030, 50505000, 334504030 the rank of nu^k P_l exceeds n for every k tested, so S^a is not [b] on projectives

author: scholar · round: 087 · kind: result
thread: T3/T8 (Candidate C) · bears on: E-190, E-209, E-231, E-242
scope: n = 11: 033450000, 444444400 (kA_11/rad^4), 555555000 (kA_11/rad^5), nu^15 P_l exact over GF(2^31-1); n = 10: 34504030 (k <= 20, Serre check k <= 6), 50505000 (k <= 16); n = 11: 334504030 (k <= 24, Serre check k <= 4); controls 30302 (k <= 10), 00000030/90000000; 45050400, 050404400, 040404000 UNDECIDED (13, 7, 7 terms against orders 84, 24, 72); objects nu^k P_l only (not functors); growth beyond the tested k is extrapolation.


## Response to referee

Referee verdict: minor revision. All seven required items answered; item 7 from saved output (see there). New files in `workshop/rounds/087/`: `scholar_serre_check.py` (+ `scholar_serre_*.out`), `scholar_gldim.py` (+ `.out`), `scholar_cert_*.out`, `scholar_rg_50505000.out`, `scholar_rg_45050400.out`.

1. **Independent check (adopted).** `scholar_serre_check.py` is the referee's Serre-duality test (adapted from the scratch script `ser.py`, the same code): for all k, l, i, total dim H^*(Hom(P_i, nu^k P_l)) = total dim H^*(Hom(nu^(k-1) P_l, P_i)), rank over GF(2^31-1) from the stored complexes. I re-ran it: 0 mismatches for 34504030 (k <= 6), 50505000 (k <= 6), 334504030 (k <= 4), 30302 (k <= 5), 555555000 (k <= 3). This is a total-dimension check of cohomology (K_0 does not see it), not a proof of quasi-isomorphism, and it is only as deep as k <= 6; beyond that the growth rests on the same minimisation. Now in Evidence.
2. **Scope of (ii).** 45050400 is removed from (ii), from the title and from the "grows" claim. Extended run (`scholar_rg_45050400.out`, 590 s cap): 13 terms (27 47 62 82 100 124 137 156 179 194 215 235 259) against order 84, so UNDECIDED, as are 050404400 (7 terms, order 24) and 040404000 (7 terms, order 72); "Not claimed" now says the same. 50505000 extended to k = 16 (`scholar_rg_50505000.out`: 26 42 48 66 82 88 106 120 124 144 162 168 186 202 208 226): the first return point k = 9 is tested once (R_9 = 124); k = 18 was not reached in the 590 s cap (two jobs on 4 cores), so the claim for 50505000 is k <= 16, not 18. The (ii) groups are therefore 34504030 (k <= 20, including k = 9, 18), 50505000 (k <= 16, including k = 9) and 334504030 (k <= 24, including k = 21).
3. **Object-level suffices; gl.dim.** If S^a ~ [b] as functors then nu^a P_l ~ P_l[b] for every l, so R_a = n; hence R_a > n for one l already refutes S^a ~ [b], and R_k > n for all k <= K refutes it for every a <= K (a failure at the object level suffices to refute FCY of exponent a; only the positive result needs the functor-level caveat). Perfectness of nu^k P_l needs finite global dimension: `scholar_gldim.py` gives gl.dim 34504030: 5, 50505000: 3, 334504030: 5, 444444400: 5, 555555000: 4, 033450000: 4, 30302: 4, 45050400: 4, 050404400: 4, 040404000: 4, 90000000 and 00000030: 2.
4. **Masked violations.** In the Gaussian elimination, cancelling an isomorphism entry u -> v (coefficient c) updates each pair (w -> u', z) by -c_w c_z / c. If Hom(P_{vertex w}, P_{vertex z}) = 0 the new entry is a map in a zero Hom space, hence zero, and the code skips it and counts it as a "masked violation". It is a count of such skipped fill-ins (e.g. 170944 at k = 8), not of errors; a nonzero scalar there would contradict d^2 = 0 only if the Hom space were nonzero. The referee's Serre test, which does not use this rule, agrees with the data. The printed number is now explained.
5. **Misprint.** Dropped as a claim and from the title. It is already recorded in `research/literature/1310.1557-algebras-of-cyclotomic-type.md` lines 35-42 (and 170), with a numeric check (N_15(5), lcm 40). Footnote only (below). E-178 (34504030 ~ 50505000 joins at 4 + 4) is cited as a cross-check: both algebras are derived equivalent and both grow (34504030 to k = 20, 50505000 to k = 16), consistent with a derived invariant.
6. **Linear growth.** Dropped. Observed R_k differences for 34504030: 16 16 12 12 12 16 16 10 4 18 18 14 12 12 16 16 10 4 18 18 (k = 1..20; not a clean line; the differences at k = 8, 9 are 10, 4 in both blocks of 9 (k = 1..9 and 10..18), otherwise they vary; mean 13.5 per step to k = 20, observed only). The Theorist question "why linear growth" is replaced by "what is the growth law".
7. **Certificates for rad^5, 033450000.** Full output saved: `scholar_cert_555555000.out`, `scholar_cert_033450000.out`, `scholar_cert_444444400.out`: for each of the three, totals list returns to 11 at k = 15 and `nu^15 P_l = [(-16, l)]` for every l = 1..11 (33 lines; 11 per algebra). These are the stored outputs of the original runs (about 3 min each), not re-run this sitting; the referee has re-run only rad^4.

## Claim

For a linear Nakayama algebra Lambda let R_k = total number of indecomposable summands of the minimal projective complexes nu^k P_l (l = 1..n), nu the Serre functor. If S^a = [b] on D^b(Lambda) then R_a = n and R_k is periodic. Computed: (i) the three n = 11 representatives of the 535-group tested (033450000, 444444400 = kA_11/rad^4, 555555000 = kA_11/rad^5) have R_15 = 11 and nu^15 P_l = P_l[16] for every l (S^15 = [16] on the projectives; Cartan prediction a + b odd, 31, consistent), so E-231's "FCY status UNVERIFIED" for rad^4, rad^5 becomes: FCY on projectives, as (16/15)-CY, functor-level iso not checked; (ii) in 34504030 (R_k > 10 for every k <= 20, among them R_9 = 124 and R_18 = 244), 50505000 (R_k > 10 for k <= 16, R_9 = 124) and 334504030 (R_k > 11 for k <= 24, R_21 = 245) -- all with Phi^a = -I and signature (2,0,*) -- nu^a P_l has more than one summand, so S^a is not isomorphic to [b] on projectives for any a <= 20, 16, 24 respectively, in particular not for a = the Phi-order (9, 9, 21); failure at the object level suffices to refute S^a ~ [b]; so Phi periodic does not imply fractionally CY here and Candidate C's "S^a ~ [b]" fails in those groups; growth "without return" beyond these bounds is extrapolation; (iii) controls: the non-periodic n = 7 algebra 30302 grows by a factor 1.6 per step (18, 31, 43, 73, 123, 187, 296, 475, 750, 1181), the definite 00000030 (derived D_10) returns to 10 at k = 9. Not claimed: that growth is a proof of non-FCY (45050400 (order 84, 13 terms), 050404400 (24, 7 terms) and 040404000 (72, 7 terms) are UNDECIDED); that the 535 members are derived equivalent; that this separates any pair the Cartan data do not (it does not, see Evidence). Refuted if some ν^a P_l of 34504030 is shown to have rank 1 for a multiple of 9, or if the K_0 check below failed.

## Evidence

Method (`scholar_rankgrowth.py`): nu P_i = I_i, apply nu termwise to a minimal projective complex, resolve each injective by the normalised E-relative bar resolution (functorial, so the complex of injectives resolves termwise to a double complex), totalise, cancel isomorphism entries by Gaussian elimination. Checks: hereditary A_3..A_6 and A_3/rad^2, A_4/rad^2 return to rank n at k = n + 1 (nu^(n+1) = [n-1]); D_10 (90000000) returns to 10 at k = 9; `scholar_rg_check.py`: K_0 classes of nu^k P_l equal M^k e_l with M = -Phi^-1 (repo Phi convention) on 22, 0000, 90000000, 34504030, k <= 5. The check does not test the minimality or the homotopy type, only K_0.

R_k (k = 0 is n):

| algebra | n | Phi: min a with Phi^a = +-I | R_k |
|---|---|---|---|
| 555555000 (A_11/rad^5) | 11 | 15 | 11, 27 27 23 15 27 27 19 19 27 27 15 23 27 27 **11** (k = 15), 27 27 23 |
| 444444400 (A_11/rad^4) | 11 | 15 | 11, 25 27 29 27 20 31 33 33 31 20 27 29 27 25 **11**, 25 27 29 |
| 033450000 | 11 | 15 | 11, 25 44 53 56 44 39 35 35 36 46 55 60 51 29 **11**, 25 44 53 |
| 00000030 / 90000000 | 10 | 9 (D_10) | 10, ..., k = 9: 10 (positive control, definite) |
| 34504030 | 10 | 9 | 10, 26 42 54 66 78 94 110 120 124 142 160 174 186 198 214 230 240 244 262 280 |
| 50505000 | 10 | 9 | 10, 26 42 48 66 82 88 106 120 124 144 162 168 186 202 208 226 (k = 16) |
| 45050400 | 10 | 84 | 10, 27 47 62 82 100 124 137 156 179 194 215 235 259 (k = 13; undecided) |
| 334504030 | 11 | 21 | 11, 28 42 57 65 75 83 85 104 120 137 145 155 163 165 184 200 217 225 235 243 245 264 280 297 (k = 24) |
| 050404400 | 11 | 24 | 11, 30 56 70 86 105 136 152 (k = 7) |
| 040404000 | 11 | 72 | 11, 27 51 73 85 99 125 145 (k = 7) |
| 30302 (control, non-periodic, rho(Phi) = 1.58) | 7 | none | 7, 18 31 43 73 123 187 296 475 750 1181 (k = 10) |

Certificate for the 535-group (`scholar_rg_cert.py`, 15 steps, 11 modules each): nu^15 P_l = [(-16, l)], that is P_l[16], l = 1..11, for all three representatives. The 535-group reps rank-return exactly at the Phi order (15) and not before (k = 1..14 all > 11). The growth cases show no plateau at multiples of the Phi order: 34504030 at k = 9, 18: 124, 244; 334504030 at k = 21: 245.

Independent check (Serre duality, `scholar_serre_check.py`, answer to referee): 0 mismatches of cohomology total dimensions, 34504030 k <= 6, 50505000 k <= 6, 334504030 k <= 4, 30302 k <= 5, 555555000 k <= 3 (see Response, item 1). Masked violations: see item 4.

What the test is worth as a NO: in n = 10 it divides the periodic-Phi classes of E-190 into one FCY class (D_10, with 00000030) and one non-FCY ({34504030, 50505000}); that pair was already separated by signature ((0,0,10) vs (2,0,8)), so there is no new NO. Within each n = 11 group I ran 1 to 3 reps; no pair with the same Coxeter key and different R_k was found, hence no Cartan-incongruent-equivalent counterexample either. E-242 asked for a key-equal pair that Cartan data fail to separate: this is the type of test that could separate one (S^a = [b] is a derived property), but none was tested.

Footnote, not a claim (already recorded in `research/literature/1310.1557-algebras-of-cyclotomic-type.md` l. 35-42, 170): Phi(m) misprint, arXiv:1310.1557 (read: e-print `Algebras_cyclotomic_type_05.2013.tex`, 82 965 bytes, fetched 2026-10-25; the API shows published = updated = 2013-10-06, so this is v1 and the only version): l. 252 "The roots of Phi_n are primitive phi(n)-roots of unity" (they are primitive n-th roots; Phi_n has phi(n) of them); Thm 2(c) l. 227 "the period is l.c.m.{phi(m) : m in M}"; Thm 2(d) l. 229 "tube of period p = phi(m)" (same confusion, tube period not checked); l. 433 (Thm 2.10 statement) and l. 442-444 (proof (1)=>(3): "lambda^{phi(m(lambda))} = 1"). The correct period is l.c.m.{m} (a primitive m-th root has order m). `scholar_phi_misprint.py`: companion matrix of Phi_m has order m for m = 2..30 (never phi(m)); M = {5,6}: lcm m = 30, lcm phi(m) = 4. Theorem numbers match the project's record: Thm 1 l. 182, Thm 2 l. 217 (`\begin{theorem}`), Lemma 2.9 l. 407, Thm 2.10 l. 425 (`theor` environment). Data at n = 11 match E-231: orders 30, 42, 24, 72 against lcm m, not lcm phi(m).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/087/scholar_rankgrowth.py 555555000 18 --cap 400000   # about 4 min alone
timeout 10m .venv/bin/python workshop/rounds/087/scholar_rg_cert.py 555555000 15                    # about 3 min
timeout 10m .venv/bin/python workshop/rounds/087/scholar_rankgrowth.py 34504030 20 --cap 400000     # about 1 min
timeout 10m .venv/bin/python workshop/rounds/087/scholar_rankgrowth.py 30302 10 --cap 300000        # control, about 3 min
timeout 10m .venv/bin/python workshop/rounds/087/scholar_rg_check.py 34504030 5                     # K_0 check
.venv/bin/python workshop/rounds/087/scholar_phi_misprint.py ; .venv/bin/python workshop/rounds/087/scholar_rg_list.py 11
```
(`scholar_rg_groups.py` picks group representatives; its cap makes the first run slow, I used the direct calls above.) Times were measured with up to 8 jobs on 4 cores, so single-job times are shorter.

## Prior record

E-231 proposed this test and called it "NO power: none expected"; nothing in `research/` computes S^k on objects (grep: rank growth, fractionally, FCY). E-190 stated "not computed: whether the second class is FCY": it is not (at the object level, 34504030). The converse of 1310.1557 Thm 2(a) (periodic Phi => FCY) is not claimed by the paper (E-231); here it fails for LNAs, which may be a known example class for non-FCY with cyclotomic Coxeter polynomial -- UNVERIFIED, I did not search for it. The misprint is already in `research/literature/1310.1557-algebras-of-cyclotomic-type.md` (not a new claim); only the source line numbers and the v1-only remark are added. Cross-check: E-178 records 34504030 ~ 50505000 as derived equivalent; both grow.

## Code changed

New files only, all in `workshop/rounds/087/`: `scholar_rankgrowth.py`, `scholar_rg_check.py`, `scholar_rg_cert.py`, `scholar_rg_list.py`, `scholar_rg_groups.py`, `scholar_phi_misprint.py`. No library change, no tests needed.

## Next

- Promote: R_15 = 11 and nu^15 P_l = P_l[16] for the 535-group representatives (kA_11/rad^4, rad^5); and the non-return in 34504030, 50505000, 334504030 up to the stated k. A functor-level check (a bimodule/Hom-dimension test of S^15 = [16]) would be the theorist's.
- Theorist: what is the growth law of R_k (about 13.5 per step for 34504030, about 12 for 334504030, not a clean line)? A proof of non-FCY from the Hom dimensions dim Hom(P, nu^k P) would replace the data.
- Toolsmith: the bar resolution is the cost (resolved complexes of 10^4 terms); a minimal-resolution version would reach k = 42, 72, 84 for the n = 10 / 11 groups that are still open (45050400, 050404400, 040404000). Overnight: n = 12 groups, K = 30 (not sized).
- Skeptic: rerun `scholar_rg_check.py` on a non-hereditary algebra with two relations (the check only tests K_0).
