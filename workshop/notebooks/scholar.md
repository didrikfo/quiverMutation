# Scholar's notebook (after round 122)

## Believe
- H-031 steps 1+2 (substitution + trivial-part removal), non-completed: the source is in Legendrian contact homology, not representation theory.
  Etnyre-Ng-Sabloff math/0101145v3: Def 3.8 (elementary/tame), Def 3.9 (stabilisation S(A) = A * E_i, d e1 = e2), Cor 3.11 (S(A) -> A homotopy eq., explicit H),
  Sec 6.3 Move III: Lemmas 6.2 (height), 6.3, 6.4, bootstrap g_i(a_i) = a_i + H(d'a_i - Phi_{i-1} d a_i) => tame iso A ~ S(A'), finitely many steps.
  Hypotheses: (T) triangular generators (digraph g -> letters of d g acyclic); (P) d a = c b + v, v free of a, b and of letters above b.
  Base Z[t,t^-1], one vertex: the transfer to K = k^n is mine (word-level; Kalck-Yang 1803.08192v1 Lemma 2.1 is a k^n analogue for vertices). Sign of H unchecked.
  Chekanov math/9709233 = PostScript only (no gs in the container); journal Invent. Math. 150 (2002) 441-483.
  Found via Manolescu-Rozenblyum 2605.08122v1 Sec 2 ("stable tame iso => q.i., [ENS, Cor 3.11]").
- scholar_ens.py (R122): (T) holds in M1, M2, M3, M3b, M4, M5, M6. ENS cancellation: M1 -> {a*, af, f} (=E-290); M4 3 pairs, 8 letters left, 3 Theta;
  M5 3 pairs, 9 left, Theta f/m, af/m, mf/m (= E-309 counts 8, 9); M2/M3/M6 no pair; M3b 2 pairs, Theta g/a, ag/a, cubic word a*.f.f.
  CORRECTION of R121/E-315: M4 is NOT "all Theta in C"; 3 Thetas stay reduced (6 thetas, 3 targets).
- So H-031 = ENS cancellation [(T)+(P), cited + k^n transfer] + Lemma A on the reduced part. Open: (T) for general Construction-2.2 parents (guess: E-284 weight).
- Completed-only (no use for kM as is): DWZ 0704.0649v4 (Def 4.1, Thm 4.6, Lemma 4.7/4.8); KY 0906.0761v3 (Thm 2.3, Lemma 2.9, 2.10, 7.4; subsection-shared counter;
  Adv. Math. 226 (2011) 2118-2168); Van den Bergh 1008.0599v3 App. A (pseudo-compact; Cor A.5.5 minimal model; Prop 7.2.1 uses the cofibrant filtration = (T)).
  Oppermann 1504.02617v1: our setting, no general lemma (Thm 7.5, Prop 7.6).
- H-024: Bastian-Holm-Ladkani 0906.3422v4 Prop 2.1 (derived eq. => P C_A P^T = C_B, any finite-dim algebras, no gldim hypothesis), Prop 2.2 (a)=>(b)=>(c)=>(d), Def 2.3.
  1012.4661v2 (J. Algebra 410 (2014)): Prop 1.3, 1.4 prop-asymmetry; Ex 2.31 ex:D15 = published Cartan-level NO (asymmetry not similar over F_3); Ex 2.21, 2.30.
  Discrepancy: 1012.4661 cites the E paper's proposition as "Prop 2.6"; counter on v4 gives 2.1 (journal numbering unchecked). Bocian-Skowronski Prop 1.5 UNVERIFIED.
  Earlier: Ladkani 1001.4765 v1 Lemma 3.5 + Aihara-Iyama 1009.3370 v3 Thm 2.27 (E-308); math/0610685 Cor 3.13, Lemma 3.14, Cor 3.15.
- Fosse 2305.06642 v1: Thm 4.2, Cor 4.3, Lemma 4.4, Thm 4.5, Cor 4.6 (rows 3,4 misprint). 2602.04571 v1 Prop 2.3; 2609.39536 v2 Thm 4.2, 7.2. Hubery-Krause 1310.1907 Cor co:complete.
- UNVERIFIED: Nikulin Cor 1.13.3 (29 secondaries now, incl. Dolgachev alg-geom/9502005 Prop 5.1 which cites 1.12.3/1.12.4/1.14.4 only); Anick 3.1/3.2; Forre 2.6; Bruning 5.1.
  Not on arXiv: Happel-Seidel, BGS, Barot-de la Pena, Nikulin 1979, Anick 1982, Forre, Bruning, CB1992, Ringel1994, Bocian-Skowronski 2005.
- Tools: curl -sL https://arxiv.org/e-print/ID into own scratch dir (gzip single tex or tar); versions/journal refs via export.arxiv.org/api/query?id_list=ID;
  search: search_query=abs:%22phrase+words%22 works, quoted multiword with %20 did not. Numbers: workshop/rounds/110/scholar_thmnumbers.py
  (python -I; [section] counters; replace \mylabel by \label first); [subsection]-shared: workshop/rounds/121/scholar_kynumbers.py.

## Did
- R106, R110, R117: E-272/273/308 era. R121: H-031 by citation negative (completed sources only), E-315.
- R122: found ENS (non-completed), tested (T)/(P) on 7 models with plants + E-309 controls (scholar_ens.py); BHL prior art for H-024; Nikulin null (+1 source).
- R122 revision (theorist, minor): cite PROOFS (Cor 3.11 proof for H; Thm 3.10 Move III proof, Sec 6.3, for bootstrap) - lemma not stated in ENS.
  H sign checked (scholar_hsign.py): left Leibniz -> (-1)^{|prefix|+1}; engine (later-letter sign) -> (-1)^{|suffix|+1}; first or last e-letter both fine.
  E-309 match = counts only (8/9 letters, 3/3 Theta). Plants = predicate power (not dg). S(M1) control added. BHL Prop 2.2/Phi need C invertible over Q.
  Asymmetry-over-F_3 (BHL) != q_S of C + C^T (E-265). M4 error was R121 scholar item 3 + R121 review item 2, not E-315's text.

## Next
- Get the R122 referee's verdict on the k^n transfer and the H sign; if asked, machine-check ENS's H on a small model (dH + Hd = tau - id on words of length <= 4).
- Ask for a PDF of Chekanov 2002 (standalone destabilisation lemma?); look in Ekholm-Honda-Kalman 1212.1519 or Ng's papers for a multi-vertex (link/composable) version over k^n.
- H-025 novelty retry (Herschend-Iyama, Lenzing-de la Pena, Kussin-Lenzing-Meltzer); prior-art line for H-026 odd n.
- Owed literature files: 0906.0761, 0704.0649, math/0101145, 0906.3422(+1012.4661), math/0702060, 2606.15595, 2604.25789, 1310.1907, 1101.3233.
- Lessons: look outside the field (symplectic DGAs solved the non-completed case); a citing paper's bibliography finds the theorem number to check;
  re-test my own earlier classifications with a new instrument (M4 was wrong); counter beats memory; separate matched from cited; YES control + planted NO.
