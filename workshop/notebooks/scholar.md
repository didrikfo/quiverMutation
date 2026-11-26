# Scholar's notebook (after round 139)

## Believe
- Λ(n,n-1) = N_n(n-1) ~der kD_n, n ≥ 4: THREE routes, none Happel-Seidel.
  (a) Brüstle 1.1+1.2 (PDF) + PD Euler form det 4, all n (E-352, rounds/137/scholar_dn.py).
  (b) Scherotzke arXiv:0906.4987v1 Thm 5.7 (LaTeX l. 1145 + PDF), proof via her 0902.0518v3 Thm 4.9, 4.15,
      and 4.15 uses BGS Thm A,B (journal, UNVERIFIED). Proof's last line typo "kA_n" for kD_n.
  (c) Happel 1987 Thm 5.12 as quoted in Hubery-Krause 1310.1907v2 Thm 7.4 (alg. closed k, simply connected)
      + explicit K0 isometry, n = 4..7 only (rounds/139/scholar_isometry.py). n ≥ 5 witness shape:
      P^-1 cols −(e1+..+e_{n-1}), e_n, e1−e_n, e2..e_{n-2}; general n unproved (asked theorist).
  With Fosse 2305.06642 Cor 4.6 (B_n ~ A_n) this is all of F-056's secondary input. LMR l. 836 only cites HS.
- BGS is NOT "cited by nothing" any more: route (b) needs it. Still not load-bearing (route (a) avoids it).
- Hubery-Krause 1310.1907v2: \begin{appendix} => numbers A.x/B.x in PDF; label co:complete = Cor A.7,
  th:subcat = Thm A.4, co:extend = A.8 (CB1992 L1), th:braid = 4.7 (CB/Ringel, HK alt proof via IS2010 = Thm 3.3),
  th:Gabriel = 4.9 (Dynkin: every indecomposable exceptional), th:Happel = 7.4. Module-level only.
  E-308 sub-steps "D^b indecomposable = shifted module" (Happel 1988) and "silting summands order into an
  exceptional sequence" stay UNVERIFIED.
- H-025 -> F-058 (R136): two-shift certificate rules out FCY; details in rounds/134, 137 notes and E-entries.
- LMR 2112.15587v2: Thm 4.4 (l. 860) conjecture-free PH classification; Thm 6.3 lead-in under Conj 6.2, 5.7;
  (2,3,10) case proved. KLM 1203.5505v1 Thm 6.11 = HS symmetry, proved.
- Primaries: Nikulin retire; Happel-Seidel replaceable (table LMR 4.4, symmetry KLM 6.11, D_n above);
  Anick load-bearing for H-023 only; BGS dependency of one route.
- Record erratum (R137, F-044 N_9(5) row) reported.
- Still UNVERIFIED: Anick 3.1/3.2; Forre 2.6; Bruning 5.1; CB1992; Ringel1994 (quoted in HK 4.7); Happel 5.12
  (quoted in HK 7.4); BGS A/B; Nikulin 1.13.3.

## Tools
- rounds/110/scholar_thmnumbers.py FILE.tex: counter-order numbering. DOES NOT handle \appendix -- check the PDF
  (curl arxiv.org/pdf/ID; pdftotext available; pypdf NOT installed in .venv).
- rounds/139/scholar_scherotzke.py NMAX: K0 shift parity Φ^{n-1} = −σ for N_n(n-1); controls A_n, E_6-8.
- rounds/139/scholar_isometry.py NMAX: backtracking isometry of Euler forms to kD_n (n=7 ~2 min, sympy).
  Its controls are weak (root counts differ); a same-root-count control is wanted.
- rounds/137/scholar_dn.py, scholar_citedeps.py; rounds/134 two-shift tools.
- arXiv: e-print is gzip of one .tex (0906.4987, 0902.0518, 1310.1907) or a tar (2112.15587); unpack in a new
  scratch dir. API: export.arxiv.org/api/query?search_query=... or id_list=ID (gives latest vN).

## Did
- R106-R134: E-272/273/308, E-315, ENS, BHL, H-027 certificate, E-339, two-shift certificate.
- R137: LMR summary confirmed; primaries table; Brüstle D_n route (E-352); KLM file drafted.
- R139: F-056 secondary input sourced (Scherotzke 5.7); K0 isometry route; wrote
  research/literature/1310.1907-categorification-non-crossing-partitions.md (v2).

## Next
- Chair items (R139 submission): F-056 label swap; E-308 "co:complete" -> A.7; README row for 1310.1907.
- Owed literature files: 0906.4987 (short), 1101.3233, math/0702060, 2606.15595, 2604.25789, LdlP math/0611532,
  Lenzing 1909.10362, math/0101145, 0906.3422; KLM 1203.5505 (draft in rounds/137/scholar.md).
- Lessons: an arXiv search by title words misses theorems buried in AR-quiver papers; search abstracts for the
  object ("path of length n-1"), not the classification. Before saying a primary is "cited by nothing",
  check the proofs of the replacement sources too (Scherotzke -> BGS). Appendix numbering: always check PDF.
