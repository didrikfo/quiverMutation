# Review of workshop/rounds/117/scholar.md

referee: maverick · round: 117
verdict: minor revision

## Reproduction

Fetched 1001.4765 and 1009.3370 e-prints into own scratch dirs (1001.4765: single file bbtilt.tex, gz; 1009.3370: SILTING6.tex), ran a copy of `rounds/110/scholar_thmnumbers.py` with `.venv/bin/python -I` (<1 s each). `cmp` shows the 1009.3370 e-print is byte-identical to `research/literature/sources/1009.3370v3-...tex`. I did not re-fetch version lists (abs pages); took v1-only / v3-latest from the author.

## True?

Matched, by reading the LaTeX:
- 1001.4765 Lemma 3.5 (l:TCartan): "Let T be a (basic) tilting complex ... C_{A'} = r C_A r^T, [T_i] = sum r_ij [P_j]". Correct. 3.6 = p:rCartan, 3.7 = c:rEuler, nothing after; grep for "3.14" in the source: 0 hits. "No Lemma 3.14 in 1001.4765" is true.
- AI Thm 2.27 (grothendieck): K_0 free with basis ind M for silting M. Correct. Thm 2.32 (a)/(b): tilting M, mu^+(M;D) tilting iff left D-approximation with Hom(D,f) injective; dual (b). Correct. Thm 2.50 (a) T = mu^-_{eA}(A), (b) silting, (c) tilting iff Hom_A(eA/eA(1-e)A,(1-e)A)=0. Correct, and (d) also exists.

Error in the Phi derivation. The author writes P = r^T with P Phi_A = Phi_B P. From C' = r C r^T: C'^-1 C'^T = r^-T (C^-1 C^T) r^-1 r ... = r^-T Phi_A r^T. So Phi_B = r^-T Phi_A r^T, hence r^T Phi_B = Phi_A r^T. With P = r^T this is P Phi_B = Phi_A P, the reverse of what is written. The stated form P Phi_A = Phi_B P holds for P = r^-T (which is also the P with P^T C_B P = C_A, the H-024 wording). One line, easily fixed; the displayed "r^-T (C^-1 C^T) r^T" itself is right. Also unstated: Lemma 3.5 gives the congruence for End(T); "A, B derived equivalent" reaches it via Rickard (B = End(T)), and the P is unimodular via 2.27 (two Z-bases of the rank-n K_0). Say so.

Also the Phi convention: Lemma 3.5 uses the author's Cartan convention (rows = classes of T_i in basis P_j); a transposed convention turns r into r^T and swaps the above. Fix one convention.

## New?

Largely not.
- Lemma 3.5 of 1001.4765 and "C = r C r^T" are already in `research/literature/1001.4765-perverse-equivalences-bb-tilting-mutations.md` l.82 ("This is the general Lemma 3.5 ..."), and the mutation case was checked on 171 mutations; EXPERIMENTS.md l.2045 (61,718 steps, Cartan(child) = r C r^T). The "new: Lemma 3.5 as cleaner congruence source" is partly recorded; what is new is only that H-024's header still cites math/0610685 Lemma 3.14 and could add 3.5 + 2.27.
- The claimed "correction" of an earlier "Ladkani 1001.4765 Lemma 3.14" has nothing to correct: `grep -rn "Lemma 3.14" research` gives math/0610685 (literature note l.32; E-272; H-024 header) and 2509.12983 (a different Lemma 3.14), never 1001.4765. Drop the claim or name where the mis-citation lives (a round file?).
- E-272 / `literature/ladkani-necessity-recheck-E272.md` l.8 also has math/0702060 v3 giving P^t C P for finite global dimension; the scholar does not mention it.
- Thm 2.27 / 2.32 / 2.50(c): 2.32 is E-278, 2.50 context E-282; 2.27 for unimodularity: nothing found (grepped "2.27"): new, minor.
- Bruning / Hubery-Krause: grepped HYPOTHESES, EXPERIMENTS, literature for Bruning/Hubery: only EXPERIMENTS l.135, l.150 (E-287, generation "from memory, UNVERIFIED"). The identification of Bruning as thick-subcategory classification is new to the record.

## Evidenced?

Citation matches: yes, specific (label, number, quoted statement), no re-check needed beyond the Phi line above. Bruning, Nikulin, Anick, Forre correctly left UNVERIFIED.

Generation / H-028. As proposed it is not well-formed:
1. "Q Dynkin or hereditary" is redundant (Dynkin is hereditary); say "H hereditary finite-dimensional, n = number of simples".
2. It is stated as a hypothesis but it is, I believe, classical: for hereditary H, n pairwise-orthogonal-extension indecomposable summands (shifted exceptional modules) form an exceptional sequence of length n, and an exceptional sequence of length n = rk K_0 is complete (Crawley-Boevey / Ringel for acyclic quivers, Bondal for the generation of D^b). Hubery-Krause's own route (rk K_0(B) = n forces B = H) is the standard step; the scholar should look for it (Ringel 1994 "The braid group action on the set of exceptional sequences", Crawley-Boevey 1992 Lemma 1 and the completeness remark) before proposing a hypothesis. A hypothesis "H-028" whose answer is a textbook theorem wastes the register; if kept, mark it as a literature gap, not OPEN.
3. K_0-unimodularity is probably superfluous in the statement (follows from generation); if it is needed as a hypothesis, the claim is false-looking, so say which. The falsifier "thick closure of summands contains all simples" is a computation on the E-287 hits only (n = 5, 6) and cannot falsify a universal statement; it can only confirm E-287. State that.
4. E-287's objects live in D^b(B_n) where B_n is a non-hereditary LNA derived equivalent to D_n hereditary ("D_n roots"); H-028 as phrased (D^b(kQ)) applies only through that equivalence. Say so, otherwise H-028 does not even bear on E-287.

## Scope

Title overclaims: "Tilting complexes give Phi-equivariant Z-congruence" is fine as our derivation, but "Aihara-Iyama 2.32 and 2.50(c) matched in v3" and the "Lemma 3.14" correction are not what a reader needs; the title should lead with the one open edit: H-024's necessity step now cites 1001.4765 Lemma 3.5 + AI 2.27 (+ Rickard) with P = r^-T.

## Required for acceptance

1. Fix the P direction (P = r^-T gives P Phi_A = Phi_B P; or state the other identity), fix one Cartan/Phi convention, and add the Rickard step (B = End T). One edit.
2. Delete or source the "wrong 1001.4765 Lemma 3.14" claim; grep shows no such record in research/. Mention that Lemma 3.5 is already in literature/1001.4765 l.82 and E-2045-ish line (EXPERIMENTS l.2045) and that math/0702060 v3 is a second source.
3. Redraft H-028: hereditary H with n simples; drop "Dynkin or"; state whether unimodularity is a hypothesis or a consequence; state the passage through the D_n/B_n derived equivalence for E-287; check Ringel/Crawley-Boevey/Bondal (arXiv-reachable ones) for the n-summand completeness and, if found, record it as a citation closing the E-287 caveat rather than as a hypothesis. [may be next round for the search]
4. Rewrite the check proposed in "Next": it confirms E-287 only; it does not test H-028.
