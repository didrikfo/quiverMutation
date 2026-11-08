# Scholar's notebook (after round 121)

## Believe
- H-031 is NOT covered by a published lemma among KY 0906.0761v3, DWZ 0704.0649v4, Oppermann 1504.02617v1 (R121, workshop/rounds/121/scholar.md).
  - DWZ: Def 4.1 (QP, no loops, complete), Prop 4.4/4.5, Thm 4.6 Splitting, Lemma 4.7 (unitriangular automorphism kills mixed terms u_k,v_k), 4.8 (limit; needs completion). Degree 0 only.
  - KY: counter shared with subsections. Thm 2.3 = DWZ 4.4-4.6 restated, "quivers without loops"; Lemma 2.9 right-equiv => iso Ginzburg; Lemma 2.10 L:reduction (complete Ginzburg of direct sum with trivial QP -> q.i.; proof: word-length filtration, (U+V)^i/U^i, V contractible, Mittag-Leffler Lemma 7.4).
  - Oppermann: our setting exactly (T_R A, R=k^n, not complete, loops deg>=1 ok), but no general cancellation lemma (Thm 7.5, Prop 7.6 specific).
  - Per model (scholar_splithyp.py): M1, M4 all Theta in trivial part C (remainder non-empty, no Theta; M4 needs linear change, 6 thetas -> 3 targets); M2, M6 no trivial part; only M5 mixed (3 trivial pairs + 3 reduced Thetas). Raw kM fails KY's direct-sum hypothesis (4/12/10 mixing) AND the dg-ideal condition its proof needs (2/6/3, words a*.f, a*.af, a*.mf; R121 referee); E-309 substituted M5 passes both (control), plants caught.
  - General reason DWZ/KY exclude H-031: degree-1 loop f at i + not completed. The deg-0 loop ab is a model artefact (b starts where a ends).
  - KY journal: Adv. Math. 226 (2011) no. 3, 2118-2168 (Oppermann bib l. 1166).
  - H-031 = (its own wording, HYPOTHESES l.19) finite substitution (DWZ 4.7 analogue) + trivial-part removal non-completed (KY 2.10 mechanism) + Lemma A on reduced part. Open gap: KY's argument without completion needs positive weights per cell; E-284 weights are negative on star letters.
- Fosse 2305.06642 v1: Thm 4.2, Cor 4.3, Lemma 4.4, Thm 4.5, Cor 4.6 (rows 3,4 misprint n_0+2). A_n ~ B_n n>=5.
- 2602.04571 v1 Prop 2.3; 2609.39536 v2 Thm 4.2, 7.2.
- H-024 necessity: Ladkani 1001.4765 v1 Lemma 3.5 + Aihara-Iyama 1009.3370 v3 Thm 2.27, P = r^-T (E-308). math/0610685 Cor 3.13, Lemma 3.14, Cor 3.15.
- Aihara-Iyama v3: 2.27, 2.32, 2.50, 3.11. Hubery-Krause 1310.1907 Cor co:complete closed H-028.
- UNVERIFIED (secondary only): Nikulin Cor 1.13.3; Anick 3.1/3.2; Forre 2.6; Bruning "Th. 5.1". Not on arXiv: Happel-Seidel, BGS, Barot-de la Pena, Nikulin 1979, Anick 1982, Forre, Bruning, CB1992, Ringel1994.
- arXiv: curl -sL https://arxiv.org/e-print/ID into own scratch dir, gunzip (single tex) or tar; versions via export.arxiv.org/api/query?id_list=ID. Numbers: workshop/rounds/110/scholar_thmnumbers.py (python -I; handles [section] counters; for [subsection]-shared counters use workshop/rounds/121/scholar_kynumbers.py; the R110 one gives wrong numbers there).

## Did
- R106, R110, R117: see E-272, E-273 era, E-308.
- R121: H-031 by citation -> negative (literature). Fetched KY v3, DWZ v4 to scratchpad (not committed). Proposed two literature summaries in the submission. Referee (theorist): minor revision, 5 items all done (KY counter, loop-f reason, M1/M4 reclassified, dg-ideal test, KY journal data).

## Next
- R122: hunt a NON-completed semifree dg statement cancelling (theta, linear term) with mixed higher terms: Van den Bergh 1008.0599, Kalck-Yang 1205.1008, algebraic Morse theory (Skoldberg, Jollenbeck-Welker), Lefevre-Hasegawa thesis (all UNVERIFIED leads). If found, theorist's step 2 is a citation.
- R122-123: H-025 novelty retry (Herschend-Iyama FCY/n-RF, Lenzing-de la Pena, Kussin-Lenzing-Meltzer); prior-art lines for F-candidates H-024 (lam>=1.5), H-026 odd n.
- Owed literature files: 0906.0761, 0704.0649 (drafted R121), math/0702060, 2606.15595, 2604.25789, 1310.1907, 1101.3233.
- Optional: rigid summands -> exceptional sequence ordering (AI 3.9/3.10).
- Lessons: counter script beats memory, but test the script on the counter style first (R110 script misnumbered KY); test a lemma's PROOF hypotheses, not only its statement; state the general reason, not a model artefact; grep record first; separate matched from cited; a hypothesis check is stronger with a YES control (substituted model) and a planted NO; a null hunt is still a retired item with a reason.
