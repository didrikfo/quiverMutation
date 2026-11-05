# Tilting complexes give Phi-equivariant Z-congruence by Ladkani 1001.4765 Lemma 3.5 + Aihara-Iyama Thm 2.27; Aihara-Iyama 2.32 and 2.50(c) matched in v3; Bruning and Anick/Forre stay UNVERIFIED

author: scholar · round: 117 · kind: backlog
thread: T3/T8, T10 · bears on: H-024, H-027, H-023 (E-265, E-267, E-272, E-273, E-278, E-282, E-287, E-289)
scope: literature only; arXiv LaTeX sources, theorem numbers by counter script (`workshop/rounds/110/scholar_thmnumbers.py`, copy used). No journal/PDF primaries. Not covered: Nikulin 1979, Anick 1982, Forre 2011, Bruning 2007 themselves.
power: n/a
backlog: item 6 settled (necessity step: matched + one-line derivation; Nikulin 1.13.3, Anick, Forre: retired as UNVERIFIED, reasons below); item 7 settled for Aihara-Iyama 2.32 / 2.50(c), generation: converted (H-entry proposed below)

## Claim

(a) The Z-congruence half of H-024's necessity step is Ladkani arXiv:1001.4765 v1 (only version) Lemma 3.5 (label l:TCartan): for a basic tilting complex T of D^b(A) with summands T_i, C_{A'} = r C_A r^T, r_ij from [T_i] = sum r_ij [P_j]. r is in GL_n(Z) because the classes of the summands of a tilting (silting) object are a Z-basis of K_0: Aihara-Iyama arXiv:1009.3370 v3 (latest) Thm 2.27. Phi-equivariance then needs no "inside a proof" citation: with Phi = C^-1 C^T (or its inverse/negative, any fixed convention built from C^-1 C^T), C' = rCr^T gives C'^-1 C'^T = r^-T (C^-1 C^T) r^T, i.e. P Phi_A = Phi_B P with P = r^T. This is our one-line algebra on a matched lemma; the statement "Phi-equivariant" is not printed in 1001.4765. The earlier citation "Ladkani 1001.4765 Lemma 3.14" is wrong: 1001.4765 has no Lemma 3.14 (section 3 ends at 3.7); Lemma 3.14 is in math/0610685 (E-272, rematched R110).
(b) Aihara-Iyama v3: Thm 2.32 ("when silting is tilting", label when silting is tilting): for tilting M, mu^+(M;D) is tilting iff every M has a left D-approximation f with Hom(D,f) injective (and dual for mu^-). Thm 2.50 (label prop:1_right_mu) (c): Okuyama-Rickard complex T = mu^-_{eA}(A) is tilting iff Hom_A(eA/eA(1-e)A, (1-e)A) = 0; (b) it is always silting. So a J != 0 step is silting, tilting only under that injectivity: E-282's "silting not tilting, by 2.32" is correctly cited.
(c) Bruning: the source is K. Bruning, Homology Homotopy Appl. 9(2) (2007) 165-176; not on arXiv, Project Euclid PDF unreachable here. Secondary statements matched: arXiv:1101.3233 v2 (proposition "Bruning", label pr:exact) and arXiv:2407.01207 v3 (proof of Prop. cites "[Bruning, Th. 5.1]"): for hereditary abelian A, A -> D^b(A) induces a bijection between exact abelian extension-closed subcategories of A and thick subcategories of D^b(A). "Th. 5.1" is therefore a secondary-source number, UNVERIFIED at the primary. Aihara-Iyama v3 does not cite Bruning (grep).
(d) Nikulin Cor. 1.13.3: only 1404.5460 (v1, v2 same line) l.336 paraphrases it ("even indefinite lattice with rank exceeding by at least 2 the minimal number of generators of d(Gamma) is determined by rank, signature, discriminant form"); no new source; UNVERIFIED stays, and no NO result uses it (E-272).
(e) Anick Thm 3.1/3.2 and Forre Th. 2.6: Forre's paper (Crelle 658 (2011) 173-192) has no arXiv version (web search and API: none); Anick 1982 not on arXiv. Secondary: Gartner 1212.2118, Efrat 2604.25789 (E-278). UNVERIFIED stays.

## Evidence

Counter script on 1009.3370v3.tex: 2.27 grothendieck, 2.28 corollary (count of summands), 2.32 when silting is tilting, 2.50 prop:1_right_mu (a mutation, b silting, c tilting iff Hom vanishes). Script on 1001.4765 main.tex: Lemma 3.5 = l:TCartan, Prop 3.6 p:rCartan, Cor 3.7 c:rEuler (finite global dimension, Euler matrix), no 3.8+; grep "congruen" in 1001.4765: 0 hits. Version check: arXiv abs pages list 1001.4765 v1; 1009.3370 v1-v3; 1404.5460 v1-v2; 1101.3233 v1-v2; 2407.01207 v1-v3 (the e-print is the latest version). Bruning searches: arXiv API (author, title words) 0 relevant hits (index thin: not absence), 5 candidate e-prints grepped, "Bruning" found only as citation in 1608.06904, 1101.3233, 2407.01207.

Generation for H-027 (E-287 "rigid, n summands, K_0-unimodular", generation UNVERIFIED). Route, steps marked: (1) over hereditary H, an indecomposable rigid object of D^b(H) is a shifted exceptional module [Happel-Ringel; UNVERIFIED]; (2) the summands order into an exceptional sequence [Aihara-Iyama Lemma 3.9/3.10 treat Hom/Ext ordering; not checked to give this exactly]; (3) K_0-unimodular n-sequence is complete in Hubery-Krause's sense (Z E = K_0; 1310.1907 v2, definitions at l.337); (4) thick(E) corresponds by Bruning to C(E) in mod H, and Hubery-Krause Thm (label th:subcat) says C(E) = mod B for a homological epimorphism A -> B with length of E = rk K_0(B); the step "rk K_0(B) = n = rk K_0(H) forces B = H" is NOT in a source read: UNVERIFIED, this is the hole. So generation is not settled, but it is not needed for "A_n ~ B_n" (Fosse Thm 4.2/Cor 4.3, E-285); it is needed only to call the E-287 chain a chain of tilting mutations.

## Reproduction

```
cp workshop/rounds/110/scholar_thmnumbers.py <scratch>/thm.py
.venv/bin/python -I thm.py research/literature/sources/1009.3370v3-aihara-iyama-silting-mutation.tex     # <1 s
curl -sL https://arxiv.org/e-print/1001.4765 (own dir, gunzip) ; .venv/bin/python -I thm.py main.tex     # <1 s each
```
Sources fetched into own directories under the scratchpad: 1001.4765, 1404.5460, 1310.1907, 1307.7573, 1608.06904, 1101.3233, 2407.01207.

## Prior record

E-272 (Ladkani math/0610685 Lemma 3.14 proof carries P; Nikulin 1.13.3 secondary), E-278 (Aihara-Iyama 2.32; Anick/Forre secondary), E-282, E-287, E-289 (generation UNVERIFIED), notebook R110. New: Lemma 3.5 of 1001.4765 as a cleaner congruence source (a recorded-then-mis-numbered "1001.4765 Lemma 3.14" is corrected), Thm 2.27 for unimodularity of r, the derivation of equivariance as algebra, Bruning identified as the thick-subcategory classification (not a generation theorem as such), and the sharpened hole in generation.

## Code changed

None. Draft literature summary: `workshop/rounds/117/scholar_literature_draft.md`.

## Next

Chair: record the Lemma 3.5 / Thm 2.27 citation in H-024 and EXPERIMENTS; fix any "1001.4765 Lemma 3.14". H-entry proposed (converted item 7): "H-028: in D^b(kQ), Q Dynkin or hereditary, a rigid object with n = rk K_0 indecomposable summands and K_0-unimodular classes is tilting (generates)"; closes with a source statement for step (4) or a short proof via Hubery-Krause (rk K_0(B)=n implies B=H) plus a computed check on the E-287 hits at n = 5, 6 (thick closure of the summands contains all simples, cheap). Theorist for the proof; toolsmith for the check. Human: Bruning and Nikulin PDFs (backlog 18).
