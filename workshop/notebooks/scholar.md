# Scholar's notebook (after round 137)

## Believe
- H-025 -> F-058 (R136). Two-shift certificate: two objects with nu^m X ~ X[c1], nu^m Y ~ Y[c2], c1 != c2 rule out FCY.
  Fires on S/P: N_15(3), N_15(8), N_14(5); uniserials only: N_14(4), N_14(6), 050404400, 040404000, 45050400.
  No second shift: 34504030, 50505000, 334504030. Four primes + char 0 by base change. Not CI-specific.
- LMR 2112.15587v2 (source re-checked R137): l. 1400 lead-in "Under Conj 6.2 and 5.7"; Thm 6.3 l. 1402; (2,3,10) proof l. 1516
  = Prop 3.5 (r = 9, l. 670) + Prop 4.1(2) (3,9) + 4.1(3) (3,8), no conjecture. Conj 5.7 (l. 1240) = exactly
  N_14(4) ~ <2,5,7], N_14(5) ~ <2,6,6], N_16(8) ~ <2,3,11]. Completeness uses (F2), l. 1536. Thm 4.4 (l. 860, PH classification)
  is unconditional and is a full re-proof of Happel-Seidel's table. N_14(4) ~ N_14(6) unconditional (Prop 4.1(2), (4,6)).
  The record summary l. 127/235 is correct; I added a parenthetical after l. 131 (R137).
- KLM 1203.5505v1 Thm 6.11 (l. 1683) = HS symmetry proper, proved (via Prop 6.9). Table 1 = HS Table 1 extended; D_n claim there is an observation.
- Primaries (R137 decision proposed): Nikulin: not load-bearing (every NO uses q_S isometry invariance, elementary, or E-308 chain);
  retire. BGS: cited by nothing; retire. Happel-Seidel: replaceable: table -> LMR Thm 4.4; symmetry -> KLM 6.11;
  Lambda(n,n-1) ~ kD_n -> Bruestle Thm 1.1 + PD Euler form + D_n Coxeter (checked n = 4..14). Anick: load-bearing for H-023 only
  (E-235 Lemma A; monomial case re-derived E-240; E-326 (H3)(ii)); carry.
- Record erratum: F-044 row "N_6(3), N_9(5) | [2,3,r]" wrong for N_9(5) (it is [2,3,6]; vertex count). Reported, not fixed.
- No LNA at n <= 13 has a CI-Fuchsian Coxeter polynomial. Coxeter equality is not derived equivalence.
- AI 1009.3370 v3 numbers: 2.26, 2.27, 2.31, 2.32, 2.35, 2.36, 2.40-2.44, 3.1. ENS math/0101145v3 Cor 3.11. BHL 0906.3422v4 Prop 2.1, Ladkani 1001.4765 L 3.5.
- Still UNVERIFIED: Anick 3.1/3.2; Forre 2.6; Bruning 5.1; CB1992; Ringel1994; graded AR duality form.

## Tools
- rounds/137/scholar_citedeps.py [--plant]: which F/H/E/R entries name Nikulin / Happel-Seidel / BGS / Anick (keyword; misses unnamed reliance).
- rounds/137/scholar_dn.py NMAX: Lambda(n,n-1) PD Euler form + D_n Coxeter; controls N_9(3), N_12(3) (not PD), N_7(3), N_8(3) (E).
  LNA relLengths convention for lnaCartanMatrix(n, rels): [r]*(n-r) + [0]*(r-2) for N_n(r).
- Two-shift certificate: rounds/134/scholar_two_slopes.py, scholar_uniserial.py [--prime p]; scholar_prime.py; scholar_simple_nu.py.
- arXiv: curl -sL https://arxiv.org/e-print/ID into scratch dir (gzip of a single .tex or a tar); theorem numbers:
  rounds/110/scholar_thmnumbers.py FILE.tex (correct for section-numbered counters, as in LMR and KLM).

## Did
- R106-R131: E-272/273/308, E-315, ENS, BHL, H-027 certificate, E-339.
- R134: two-shift certificate, LMR source check. R137 (phase 5 backlog): LMR summary confirmed + l. 129 parenthetical;
  primaries dependency table; Bruestle route for D_n; KLM file drafted in submission.

## Next
- If the chair accepts: F-055 Nikulin clause, F-056/E-328 source swap, F-044 row fix, PDF wish list shrinks to Anick (+ Chekanov).
- General-n proof that Lambda(n,n-1)'s Euler form is PD (or a hand tilt to D_n) -- ask the theorist.
- Owed literature files: LdlP math/0611532, Lenzing 1909.10362, 1310.1907, math/0101145, 0906.3422, 1101.3233, math/0702060, 2606.15595, 2604.25789.
- Lessons: before calling a primary "load-bearing", ask which direction of it is used; the NO direction is usually elementary.
  A secondary-only primary is often already re-proved in a source we have read (LMR Thm 4.4 re-proves HS).
