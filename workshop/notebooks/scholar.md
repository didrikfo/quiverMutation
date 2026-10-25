# Scholar's notebook (after round 087)

## Believe
- Candidate C = S^a(A) ~ A[b] (fractionally CY). 1310.1557 Thm 2(a) FCY => Phi periodic; the converse is NOT a theorem and (R087, rank growth) fails for LNAs: periodic Phi with signature (2,0,*) (n = 10: 34504030, 50505000, 45050400; n = 11: 334504030) has R_k = total rank of nu^k P_l growing ~ linearly with no return; the n = 11 535-group reps (033450000, kA_11/rad^4, rad^5) return to 11 at k = 15 with nu^15 P_l = P_l[16] (projective level only). Growth is data (k <= 24), not proof. 24- and 72-groups: only k <= 7 run.
- Method works: bar-resolve nu P_i, Gaussian-minimise (workshop/rounds/087/scholar_rankgrowth.py); checks: A_n, A_n/rad^2 period n+1, D_10 period 9, K_0 = -Phi^-1 powers.
- 1310.1557 is v1 only (API published = updated 2013-10-06). Misprint phi(m) for m: l. 227 (Thm 2(c)), 229, 252, 433, 442-444. Thm 1 l. 182, Thm 2 l. 217, Lemma 2.9 l. 407, Thm 2.10 l. 425 (matched R087).
- 0911.5137 (LaTeX checked R082): Thm A/B/C l. 272/363/526; Cor 1.2 l. 439, 1.3 l. 488, 1.9 l. 671 (silent at prime n = 11), 1.13 l. 967, Remark 1.14 l. 972.
- Older: AI Thm 2.31, 2.32(b); Oppermann 1504.02617 v1 Thm 1.1 / Prop 7.6 (E-227); T5 gap = repo child = H_0(kM) (E-217).
- R085 (T5): Lemma A (Anick, monomial quadratic no-overlap => strongly free) has no reachable theorem number; nearest Loefwall 1711.03644v2 Sect. 3. F-047's 3304000 is not almost separate (named via lnaMoves mates + Fosse thm:QuipuToAn).
- arXiv e-print works from the cloud (curl -sL https://arxiv.org/e-print/ID).

## Did
- R001..R085 audits and summaries. R087: rank-growth test, scripts scholar_rankgrowth/rg_check/rg_cert/rg_list/rg_groups/phi_misprint in rounds/087.

## Next
- Reach k = 42/72/84 for 045050400-type groups (faster resolution) and n = 12 groups; look for a key-equal pair with different growth (the only way this test beats Cartan data).
- Literature: is "cyclotomic Coxeter polynomial but not FCY" known for Nakayama (Herschend-Iyama, Lenzing-de la Pena)? Not searched.
- Anick 1982/86 theorem number for Lemma A; Fosse's quipuForAlmostSeparateLNA vs the theorem statement.
- Lessons: grep summaries first; remarks/props/cors share the lemma counter; never pkill -f a pattern that matches the own shell command; run jobs one per core.
