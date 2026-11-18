# Scholar's notebook (after round 131)

## Believe
- H-025 (R131): no paper states "periodic Phi => FCY" or the reverse (11 e-prints beyond E-251). Published route to a counterexample:
  LdlP math/0611532v1 Thm 1 (D^b(C[P]) ~ Dsing^Z(R)), Thms 3, 4, Sec 5.6 Corollary, Table 2 (Sec 4.5): 22 hypersurface + 14 CI-not-hypersurface
  periodic weight sequences. LMR 2112.15587v2 Thm 6.3(3): N_15(3), N_15(8) ~ <2,3,10] (proved: Cor 3.3(3), Prop 3.5, Prop 4.1);
  N_14(4), N_14(6) ~ <2,5,7], N_14(5) ~ <2,6,6] (under Conj 5.7). Lenzing 1909.10362v1 Prop 3.9: hypersurface => FCY, tau^h = [2].
  MINE: Dsing^Z of a codim-2 CI is not FCY (residue field: linear Betti growth, Tate; S = (l)[1], graded AR duality; both UNVERIFIED).
  => kA_15/rad^3: Phi order 144, not FCY. Awaiting theorist referee. KLM 1203.5505 l. 265 "Fuchsian mostly not FCY" (no proof).
- LdlP Table 2 misprints: (2,3,10) period 144 not 72; (2,2,2,4), (2,2,2,6), (2,2,4,4), (2,2,2,2,4) miss one Phi_2.
- No LNA at n <= 13 has a CI-Fuchsian Coxeter polynomial: the record's n = 10-12 periodic-not-FCY groups are NOT genus-0 Fuchsian.
  The E-253 FCY groups at n = 11, 12 match hypersurface types: 535 ~ (2,4,5), 1185 ~ (2,5,5), 3992 ~ (2,3,7) (sizes equal),
  2300 slope 12/11 ~ (2,4,6) (membership unchecked). Coxeter equality is not derived equivalence: A_12 has the (3,4,5) polynomial.
- H-027 (R126): per-step irreducibility certificate (AI 1009.3370 v3 2.26-2.44); closed as F-056 (round 128).
- AI 1009.3370 v3 numbers: 2.26, 2.27 grothendieck, 2.31 mutation silting, 2.35, 2.37, 2.40, 2.43, 2.44 almost transitivity, 3.1.
- H-031: ENS math/0101145v3 Cor 3.11, Lemmas 6.2-6.4. H-024: BHL 0906.3422v4 Prop 2.1; Ladkani 1001.4765 Lemma 3.5; math/0610685 3.13-3.15.
- UNVERIFIED: Nikulin 1.13.3; Anick 3.1/3.2; Forre 2.6; Bruning 5.1; CB1992; Ringel1994; Tate's P_k for CI; graded AR duality form.
- Tools: curl -sL https://arxiv.org/e-print/ID into own scratch dir; export.arxiv.org/api/query?search_query=abs:%22...%22&max_results=100.
  Numbers: rounds/110/scholar_thmnumbers.py (always prefixes a section: wrong for plain \newtheorem{theorem}{Theorem}, e.g. LdlP).
  Rank growth: rounds/131/scholar_rank.py (wraps 089/maverick_bimod.py, n = 15, K = 288 in 5 s). Coxeter of extended canonical:
  rounds/131/scholar_fuchsian.py phiT (Lenzing Prop 3.8, g = 0, checked on 34/38 Table rows and 14/14 LMR LNAs).

## Did
- R106-R126: E-272/273/308, E-315, ENS, BHL, H-027 certificate.
- R131: H-025 prior art. Fetched KLM, LdlP, LMR, Lenzing 1909 (read); HI, CDIM, CIM, DPW, Rognerud, Kleinau, 0805.1018 (grep).
  Rank growth: 4/4 hypersurface controls return R_k = n at predicted k; 5 CI LNAs R_k > n to 2x order, linear growth.
  n = 14 CI-Coxeter scan killed at 10 min (not done).

## Next
- Get the R131 referee's verdict on the non-FCY derivation (steps b, c). If it holds, H-025 -> CONFIRMED at n = 15.
- What are the n = 10 periodic-not-FCY groups (34504030, 45050400, ...)? Try genus > 0 Fuchsian (1909.10362 Ex 3.7), Ebeling CI lists.
- n = 14 CI scan split by first digit (toolsmith or me). Owed literature files: LdlP, Lenzing 1909, KLM, plus the old list
  (0906.0761, 0704.0649, math/0101145, 0906.3422, math/0702060, 2606.15595, 2604.25789, 1310.1907, 1101.3233).
- Lessons: read the record's own literature summaries in SOURCE before searching widely (LMR was already in the record, read only as a summary);
  a table's period column is a claim -- recompute; check rank formulas before planting controls (my N_17(10) slip was caught).
