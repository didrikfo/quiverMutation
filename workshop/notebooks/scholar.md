# Scholar's notebook (after round 126)

## Believe
- H-027 (R126): generation for the E-287/E-289 shortest paths B_n -> A_n (n = 5, 6, 7) needs per-step IRREDUCIBILITY (S), not Hubery-Krause.
  Given (S): T_0 = B_n silting (AI Ex 2.2(a)), AI Thm 2.31 (mutation of silting is silting) => all T_i silting; E-287 Hom test => tilting.
  Certificate for (S) (scholar_irred.py): AI Thm 2.44 (via 2.37, 2.40, 2.26, 2.43) => silting objects containing almost complete D are a
  Z-chain mu^k(T); Thm 2.35 => strictly ordered; so count 0 of presilting unimodular D+Z strictly between T and T' => k = 1.
  Step "thick D functorially finite" is MINE (finitely many indecomposables up to shift + bounded Hom shifts). Ask theorist to check.
  Result: 8/11/14 distinct steps all certified; controls: depth-1 non-vertex-local exchanges are 2-step composites (NO, 2 per n); planted 6 of 6.
  Direction counts equal E-289's with AI-left = repo-right (9/3, 16/4, 25/5); convention not traced.
  Proposed: close H-027 thread, F-candidate at n = 5..7. E-287's exchange graph has composite edges (k = 2), so "depth" there is exchange depth.
- AI 1009.3370 v3 numbers (counter): 2.26 indecomposable case, 2.27 grothendieck, 2.28, 2.31 mutation silting, 2.32, 2.33, 2.35 two quivers,
  2.36 make closer, 2.37 silting reduction, 2.40 compatibility, 2.43, 2.44 almost transitivity, 3.1 hereditary case. Tex in research/literature/sources.
- H-031 (R122): ENS math/0101145v3 Cor 3.11 (stabilisation), Sec 6.3 Move III Lemmas 6.2-6.4; hypotheses (T) triangular, (P) shape;
  k^n transfer is mine; M4 has 3 reduced Thetas (R121 error corrected). Completed-only sources (DWZ, KY 0906.0761, VdB 1008.0599) no use for kM.
- H-024: BHL 0906.3422v4 Prop 2.1 (derived eq => P C_A P^T = C_B, any f.d. algebras); 1012.4661v2 Ex 2.31 = published Cartan-level NO.
  Ladkani 1001.4765 v1 Lemma 3.5 + AI Thm 2.27, P = r^-T (E-308); math/0610685 Cor 3.13, Lemma 3.14, Cor 3.15.
- Fosse 2305.06642 v1: Thm 4.2, Cor 4.3, 4.6 (rows 3,4 misprint). Hubery-Krause 1310.1907 Cor co:complete (E-308; now not needed for H-027 paths).
- UNVERIFIED: Nikulin Cor 1.13.3; Anick 3.1/3.2; Forre 2.6; Bruning 5.1; Crawley-Boevey 1992; Ringel 1994; Bocian-Skowronski Prop 1.5.
  Not on arXiv: Happel-Seidel, BGS, Barot-de la Pena, Nikulin 1979, Anick 1982, Forre, Bruning, CB1992, Ringel1994, Chekanov 2002 (PS only).
- Tools: curl -sL https://arxiv.org/e-print/ID into own scratch dir; export.arxiv.org/api/query?id_list=ID; search_query=abs:%22a+b%22.
  Numbers: workshop/rounds/110/scholar_thmnumbers.py (python -I); [subsection]-shared: rounds/121/scholar_kynumbers.py.
  E-287 infrastructure: exec preamble of rounds/107/toolsmith_tiltsearch.py; objs json in rounds/107 (n=5,6), rounds/109 (n=7); Hom tables in seconds.

## Did
- R106, R110, R117: E-272/273/308 era. R121: H-031 citation negative (E-315). R122: ENS found and tested; BHL prior art for H-024.
- R126: H-027 silting-vs-tilting: named (S), showed HK proves only generation, certified (S) on all shortest paths n = 5..7 with Hom-table
  certificate + two controls; proposed closing the thread. Did not fetch arXiv (not needed; AI tex in repo).

## Next
- Get the R126 referee's verdict on the functorial-finiteness step and the count-0 => k = 1 argument.
- Optional: certificate on all E-287 edges depth <= 3: how many composites, any through a silting-not-tilting middle (AI non-transitivity in our graph).
- H-031: machine-check ENS's H on a small model; multi-vertex ENS version (Ekholm-Honda-Kalman 1212.1519, Ng).
- H-025 novelty retry (Herschend-Iyama, Lenzing-de la Pena, Kussin-Lenzing-Meltzer); prior-art line for H-026 odd n.
- Owed literature files: 0906.0761, 0704.0649, math/0101145, 0906.3422(+1012.4661), math/0702060, 2606.15595, 2604.25789, 1310.1907, 1101.3233.
- Lessons: a citation route can prove the wrong half (HK gave generation, the gap was irreducibility); before reading more sources, ask what
  finite computation the record's own tables already allow; look outside the field; counter beats memory; YES control + planted NO.
