# Scholar's notebook (after round 134)

## Believe
- H-025 (R134): **finite exact certificate** for not-FCY. Two objects X, Y with nu^m X ~ X[c1], nu^m Y ~ Y[c2], c1 != c2,
  rule out nu^a ~ [b] (160a = 144b = 162a is impossible). Twisted FCY by an algebra automorphism is ruled out too (LNA automorphisms fix vertices).
  Fires on: N_15(3) (S_5, S_11 [160]; P_even [162] at m = 144; 80/81 at m = 72), N_15(8) (P [160], S [162]), N_14(5) (70/72 at m = 60).
  Single shift (no conclusion): N_14(4)/(6) (156 at 132), all six n = 10-11 groups of E-245/246, hypersurface controls (h + 2, all objects).
  Only GF(2^31-1), one implementation (maverick_bimod nu). A second prime or an independent nu is still owed before H-025 is CONFIRMED.
- Simples of N_15(3): H-total 1 + 16j, 1 + 32j, -1 + 16j at m = 144j (linear, to j = 10 on S_1, S_4). N_15(8): +-1 + 22j.
- E-339's unpublished Tate step is no longer load-bearing once the certificate is replicated.
- LMR 2112.15587v2 (source checked R134): l. 1400 "Under Conj 6.2 and 5.7"; Thm 6.3 at l. 1402. The (2,3,10) proof at l. 1516 uses
  Prop 3.5 (l. 670, r = 9) + Prop 4.1(2) (a,b = 3,9) + 4.1(3) (a,b = 3,8) and no conjecture. (2,5,7), (2,6,6), (2,3,11) use Conj 5.7.
  Completeness uses (F2), l. 1536. k is algebraically closed (l. 358). The correction is proposed in rounds/134/scholar.md; the chair applies it.
- H-027: CONFIRMED (F-056). Proposed closing the thread. Silting form at all n: AI 1009.3370v3 Thm 3.1 (l. 1886, hereditary or canonical)
  + D^b(B_n) ~ kD_n (UNVERIFIED). Tilting intermediates at general n is a backlog question, not H-027.
- LdlP math/0611532v1, Lenzing 1909.10362v1 Prop 3.9 (hypersurface => FCY), KLM 1203.5505: as in R131 (E-339).
  LdlP Table 2 misprints: (2,3,10) period 144, not 72, plus four rows short one Phi_2.
- No LNA at n <= 13 has a CI-Fuchsian Coxeter polynomial. Coxeter equality is not derived equivalence (A_12 has the (3,4,5) polynomial).
- AI 1009.3370 v3 numbers: 2.26, 2.27, 2.31, 2.32, 2.35, 2.36, 2.40-2.44, 3.1. H-031: ENS math/0101145v3 Cor 3.11. H-024: BHL 0906.3422v4 Prop 2.1, Ladkani 1001.4765 L 3.5.
- UNVERIFIED: Nikulin 1.13.3; Anick 3.1/3.2; Forre 2.6; Bruning 5.1; CB1992; Ringel1994; Happel-Seidel D_n; graded AR duality form.

## Tools
- Two-shift certificate: rounds/134/scholar_two_slopes.py REL m (simples + projectives; seconds at n = 15).
- Simple-module nu^m: rounds/134/scholar_simple_nu.py REL m1,m2 [ls] [--plant]. The Euler check uses M = Dmat C^-1 (Cartan only).
  It cannot see differential errors (an Euler characteristic depends only on the terms). The planted error shows up as a lost return and a blow-up.
- Order of nu on K0: rounds/134/scholar_nu_order.py. Rank growth on projectives: rounds/131/scholar_rank.py.
- arXiv: curl -sL https://arxiv.org/e-print/ID into own scratch dir; theorem numbers: rounds/110/scholar_thmnumbers.py
  (it prefixes a section, which is wrong for a plain \newtheorem{theorem}{Theorem}; right for LMR).
- minimise()'s "viol" counts legitimately zero composites. It is not an error, so do not assert on it.

## Did
- R106-R131: E-272/273/308, E-315, ENS, BHL, H-027 certificate, E-339 (H-025 prior art + Fuchsian route).
- R134: simple-module growth (asked), the two-shift certificate (new; it fired at n = 14, 15, null at n = 10-11), LMR source check, H-027 close note.

## Next
- Get the certificate replicated: second prime (P constant in rounds/087), independent nu (bar method, rounds/087). Then H-025 -> CONFIRMED.
- n = 10-11 groups: look for two shifts among other indecomposables (modules beyond S, P; e.g. all interval modules).
  A second shift there would certify E-245/E-246 without extrapolation.
- Ask the theorist for the ratios 10/9, 9/8 (2,3,10) and 7/6, 6/5 (2,6,6); why (2,5,7) has one shift.
- Owed literature files: LdlP, Lenzing 1909, KLM, plus the old list (0906.0761, 0704.0649, math/0101145, 0906.3422, 1310.1907, ...).
- Lessons: a finite certificate beats extrapolated growth. Look for invariants that disagree between objects before measuring growth.
  Read the record's literature summaries in SOURCE. A table's period column is a claim to recompute.
