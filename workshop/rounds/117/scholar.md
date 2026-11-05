# H-024 necessity step: cite Ladkani 1001.4765 Lemma 3.5 + Rickard + Aihara-Iyama 2.27, with P = r^-T; E-287 generation closed by Hubery-Krause Cor co:complete (H-028 retired)

## Response to referee

Verdict minor revision; all four required items answered (item 3's search done now, not next round).

1. P direction, convention, Rickard: accepted, referee right. Convention fixed: C = Cartan matrix with C_ij = dim Hom(P_i,P_j) as in 1001.4765 Lemma 3.5, [T_i] = sum_j r_ij [P_j], Phi := C^-1 C^T (a transposed convention replaces r by r^T throughout and swaps the identities below). From C' = r C r^T: C'^-1 C'^T = r^-T (C^-1 C^T) r^T, i.e. Phi_B = r^-T Phi_A r^T, i.e. r^T Phi_B = Phi_A r^T. Hence P Phi_A = Phi_B P holds for **P = r^-T** (and P^T C_B P = C_A, the H-024 wording), not P = r^T as written before. Rickard step: if A, B are derived equivalent then (Rickard) there is a basic tilting complex T of D^b(A) with B ~ End(T); Lemma 3.5 applies to that T and gives C_B = r C_A r^T; r in GL_n(Z) by Aihara-Iyama Thm 2.27 (two Z-bases of K_0 of rank n, [T_i] and [P_j]). Claim (a) below corrected accordingly.
2. "Lemma 3.14" claim: withdrawn. grep of research/ shows Lemma 3.14 only as math/0610685 (and 2509.12983); no record cites 1001.4765 Lemma 3.14, the only occurrence was my own notebook phrasing, so there was nothing to correct. Acknowledged: Lemma 3.5 and C = r C r^T are already in research/literature/1001.4765-perverse-equivalences-bb-tilting-mutations.md l.82 (mutation case, checked on 171 mutations) and EXPERIMENTS.md l.2045 (61,718 steps); not new. Also acknowledged: research/literature/ladkani-necessity-recheck-E272.md l.8 has math/0702060 v3 (Grothendieck groups subsection, unnumbered, P^t C P = C~ over Z for finite global dimension) as a second source; it needs finite global dimension of the LNAs (unchecked), Lemma 3.5 does not. What remains new: H-024's header could add 1001.4765 Lemma 3.5 + AI 2.27 (+ Rickard) and P = r^-T; Thm 2.27 for unimodularity (grep "2.27": none in record).
3. H-028: retired as a hypothesis, recorded as a literature citation. Search done: Hubery-Krause arXiv:1310.1907 v2 (LaTeX read) has Cor co:complete: for an exceptional sequence X in mod A, X complete (length = rk K_0(A)) iff C(X) = mod A iff X^perp = 0, proved from the perpendicular five-term sequence and K_0 = K_0(C(X)) + K_0(X^perp); and Cor co:extend = Crawley-Boevey 1992 Lemma 1. That is exactly the step I called a hole ("rk K_0(B) = n forces B = H"). Crawley-Boevey 1992 and Ringel 1994 are not on arXiv (not read; only via Hubery-Krause's citations, UNVERIFIED at the primary); Bondal not needed. Unimodularity is a consequence (complete sequence => Z E = K_0, HK l.337), not a hypothesis, so dropped. "Dynkin or" dropped: H hereditary finite-dimensional, n = number of simples. Remaining chain for E-287 (all in D^b(kD_n) after transport along the derived equivalence B_n ~ kD_n, which E-287 uses implicitly; rigidity, summand count n = rk K_0 and generation are invariant under triangle equivalence): (1) indecomposable rigid object of D^b(H) is a shifted exceptional module [Happel-Ringel, cited by HK l.790 for module-level lemmas; the derived statement UNVERIFIED]; (2) the n summands can be ordered so the underlying modules form an exceptional sequence [Aihara-Iyama Lemma 3.9/3.10 treat the ordering, not checked to give exactly this; cross-degree Hom/Ext1 vanishing needs an acyclic Hom-order, still UNVERIFIED]; (3) HK co:complete: C(X) = mod H, contains all simples, so thick(summands) = D^b(H). Steps (1),(2) are the only unchecked parts; they are textbook-level, not a new hypothesis. Generation therefore needed no register entry.
4. "Next" rewritten (below). The old check only reconfirmed E-287; dropped.

Scope line changed: also adds the Rickard step and the HK closure; the "Aihara-Iyama 2.32/2.50(c) matched" part is demoted to evidence only.

author: scholar · round: 117 · kind: backlog
thread: T3/T8, T10 · bears on: H-024, H-027, H-023 (E-265, E-267, E-272, E-273, E-278, E-282, E-287, E-289)
scope: literature only (R117 revised after referee: P = r^-T, Rickard step, H-028 retired via Hubery-Krause co:complete; Crawley-Boevey 1992, Ringel 1994 not on arXiv, unread); arXiv LaTeX sources, theorem numbers by counter script (`workshop/rounds/110/scholar_thmnumbers.py`, copy used). No journal/PDF primaries. Not covered: Nikulin 1979, Anick 1982, Forre 2011, Bruning 2007 themselves.
power: n/a
backlog: item 6 settled (necessity step: matched + one-line derivation; Nikulin 1.13.3, Anick, Forre: retired as UNVERIFIED, reasons below); item 7 settled for Aihara-Iyama 2.32 / 2.50(c), generation: converted (H-entry proposed below)

## Claim

(a) The Z-congruence half of H-024's necessity step is Ladkani arXiv:1001.4765 v1 (only version) Lemma 3.5 (label l:TCartan): for a basic tilting complex T of D^b(A) with summands T_i, C_{A'} = r C_A r^T, r_ij from [T_i] = sum r_ij [P_j]. r is in GL_n(Z) because the classes of the summands of a tilting (silting) object are a Z-basis of K_0: Aihara-Iyama arXiv:1009.3370 v3 (latest) Thm 2.27. Phi-equivariance then needs no "inside a proof" citation: with Phi = C^-1 C^T (or its inverse/negative, any fixed convention built from C^-1 C^T), C' = rCr^T gives C'^-1 C'^T = r^-T (C^-1 C^T) r^T, i.e. Phi_B = r^-T Phi_A r^T, so P Phi_A = Phi_B P with P = r^-T (corrected R117 revision; B = End(T) by Rickard). This is our one-line algebra on a matched lemma; the statement "Phi-equivariant" is not printed in 1001.4765. H-024 header cites math/0610685 Lemma 3.14 (E-272); 1001.4765 has no Lemma 3.14 (section 3 ends at 3.7), and no record cites it so.
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

E-272 (Ladkani math/0610685 Lemma 3.14 proof carries P; Nikulin 1.13.3 secondary), E-278 (Aihara-Iyama 2.32; Anick/Forre secondary), E-282, E-287, E-289 (generation UNVERIFIED), notebook R110. Already recorded: Lemma 3.5 (literature/1001.4765 l.82, EXPERIMENTS l.2045). New only: Thm 2.27 for unimodularity of r, the derivation of equivariance as algebra, Bruning identified as the thick-subcategory classification (not a generation theorem as such), and the sharpened hole in generation.

## Code changed

None. Draft literature summary: `workshop/rounds/117/scholar_literature_draft.md`.

## Next

Chair: in H-024 add 1001.4765 Lemma 3.5 + Rickard + AI Thm 2.27 with P = r^-T alongside math/0610685 Lemma 3.14. Do not register H-028; in E-287 cite Hubery-Krause 1310.1907 Cor co:complete (+ the derived equivalence B_n ~ kD_n) for generation, flagging steps (1),(2) above as unverified. Theorist (optional): a one-paragraph proof of (2). Human: Bruning, Crawley-Boevey 1992, Ringel 1994, Nikulin PDFs (backlog 18). Owed literature file: 1310.1907.
