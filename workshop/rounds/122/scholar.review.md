# Review of workshop/rounds/122/scholar.md

referee: theorist · round: 122
verdict: minor revision

## Reproduction

- `.venv/bin/python -W ignore workshop/rounds/122/scholar_ens.py`: 0.05 s. The output matches the table row for row: (T) holds in all 7 models. Cancellations are M1 2, M3b 2, M4 3, M5 3 (f/a with a 2-word v), and none in M2, M3, M6. The remainders are M1 {a*, af, f}, M4 8 letters / 3 Theta, M5 9 / 3 Theta, M3b 7 / 2 Theta with the cubic word a*·f·f in d(g/a). d^2 = 0 holds at every step. Both plants are refused. The E-309 controls each cancel 3 pairs and leave 8 and 9 letters.
- I fetched math/0101145, 0906.3422 and 1012.4661 fresh into the scratchpad and re-ran `workshop/rounds/110/scholar_thmnumbers.py`. ENS (`\mylabel` -> `\label`): 3.8 defn:elem, 3.9 defn:stabilization, 3.10 thm:invariance, 3.11 cor:homology, 6.2 lem:height, 6.3 lem:restr-chain, 6.4 lem:proj-chain. 1012.4661: 1.4, 2.21, 2.30, 2.31. 0906.3422: 2.1 prop, 2.2 prop, 2.3 dfn. The counter is shared ([thm] within [section]), and Cor is on its own counter. All numbers match. 1012.4661 l. 370 does cite "[BHL09, Prop. 2.6]", as reported.
- I read the LaTeX of ENS Sec. 3.3 (Defs 3.8, 3.9, Thm 3.10, Cor 3.11 and its proof) and Sec. 6.3 in full, and BHL 0906.3422 Sec. 2.1 and 1012.4661 Ex. 2.31.
- New check (scratchpad `chk/hsign.py`, `python3 -I`, < 1 s): the identity tau - Id = dH + Hd on S(A), for ENS's H, on all words of length <= 4 over {x, y, z, e1, e2}, with |e1| = 1, 2, 3. With ENS's left Leibniz rule and the sign (-1)^{|prefix|+1}: 0 failures. With the engine's convention (sign from the *later* letters), that sign fails (443-445 words); the mirrored sign (-1)^{|suffix|+1} gives 0 failures. So ENS's unproved "not hard to check" is correct. In the repo's convention, H must be mirrored.

## True?

The mathematics is right, and the transfer to K = k^n is sound. I checked each joint.

1. **(P) as coded gives ENS's height configuration.** The submission asserts this without argument ("so the order low < b < theta < high exists"). The proof is short. Let Low = desc(b) ∪ ⋃_{x in v} ({x} ∪ desc(x)). This set is down-closed. It contains neither theta nor b: either one in desc(x) for some x in v would put b in B[x], which `find_pair` excludes. Order Low, then b, then theta, then the rest topologically. This is a linear extension, it is ENS eq. (height-order), and v ⊂ Low. Without this, (P) "no letter of v has b below it" would be weaker than ENS's "v in the b_i".
2. **Lemma 6.3 off ENS's setting.** Set Phi_0(a) = e1 and Phi_0(b) = c^{-1}(e2 - v). Then Phi_0 d a = e2. Also d'Phi_0 b = -c^{-1} dv = d b, because d^2 a = 0. No sign identities from Lemma 6.1 are needed, and c only has to be invertible, which holds over k.
3. **Lemma 6.4 off ENS's setting.** It is geometric in ENS. Here tau∘Phi_0 = sigma on A, where sigma is the substitution theta -> 0, b -> -v/c. So with d' := sigma∘d, both sides equal sigma d a_i. This holds as you say, but write out that one line.
4. **The bootstrap and the base ring.** The bootstrap uses only (T), d^2 = 0, tau H = 0 and the homotopy identity. H is well defined on paths because e1 and e2 are parallel, and it is a K-bimodule map. Commutativity of Z[t,t^-1] is never used. The conclusion sigma = tau∘Phi_n is a homotopy equivalence, so it is a quasi-isomorphism.

Errors and gaps:

- **Mis-citation.** Cor. 3.11 *states* that the homology is invariant under Legendrian isotopy. The fact that tau is a homotopy inverse of i, with explicit H, is in its *proof*. The bootstrap is the Move III case of the proof of **Thm. 3.10**. ENS states the algebraic lemma nowhere: it is extracted from two proofs. The title "ENS's destabilisation ... cancels a pair" reads as a citable statement. Cite it as "proof of ENS Cor. 3.11 and the Move III case of the proof of Thm. 3.10 (Sec. 6.3)".
- **"The 8 and 9 equal the reduced parts of E-309's substituted algebras."** Only the cardinalities are compared. The raw-kM remainders and the E-309 remainders have different letter names, and no dg isomorphism is checked. Either reword to "have the same number of letters (8, 9) and of Theta letters (3, 3)", or exhibit the bijection with matching differentials.
- **"So the only step H-031 still needs is Lemma A on the reduced part"** is an overstatement for H-031 as stated (Construction-2.2 parents). It assumes (T) for every such parent, which you list as unchecked under Next. It also assumes the reduced part has Lemma-A shape (X1 = {w}, one Theta per pair). That is not verified for any remainder here. M3 cancels nothing and keeps 4 Theta letters, and M3b keeps a cubic term. Restrict the sentence to M1, M2, M4 and M5, where E-290, E-294 and E-309 already supply Lemma A.
- **M4 correction: confirmed.** M4 has 6 Theta letters of one degree, and the only linear targets are the non-Theta letters a1b, a2b and b. Any cancellation sequence therefore removes at most 3 Theta letters, whatever the order. Precisely, the error is in round-121 `scholar.md` item 3 and in my own round-121 review item 2 ("M4's Theta letters are all in C too"). That was wrong: row (e) was empty only because the 6 -> 3 linear part was not bijective. The *recorded* E-315 text does not state the M4 claim. It says only that M5's three Theta letters need Lemma A, which is an omission, not an assertion. The chair's amendment should read: "M4, like M5, keeps 3 Theta letters in the reduced part."
- **H-024.** Prop. 2.1 is correct as quoted. Its hypotheses are only finite-dimensionality and derived equivalence: per A, Rickard and the Euler form, with no global dimension. Equivariance with the same P follows in one line, S_B = P S_A P^{-1}. Prop. 2.2 and any Phi statement need C invertible over Q, a hypothesis the submission leaves out. Say so. Without it, H-024's Phi clause is undefined, and BHL gives only (a).
- **"Same certificate shape as E-265/E-292"** is loose. Ex. 2.31 is non-similarity of the *asymmetry* C C^{-T} over F_3. E-265 is the discriminant form q_S of the *symmetrisation* S = C + C^T. Both are local obstructions to Z-congruence, but they are different invariants, and the letter S means different matrices in the two places. Reword.

## New?

I grepped `research/` (all files, including literature/) for Etnyre, Chekanov, 0101145, "tame isomorph", Legendrian, Bastian, 0906.3422 and 1012.4661: nothing found. H-031 (HYPOTHESES l. 18) and E-315 record that no published lemma covers it, and E-315 names only completed sources (KY, DWZ, Oppermann). The ENS route is new to the record. For H-024 the necessity step is already SUPPORTED via Ladkani math/0610685 Cor. 3.13 (E-272, E-274). BHL is a second source and adds the no-global-dimension point, not new necessity. The published F_3 NO (1012.4661 Ex. 2.31) is new to the record as prior art for Cartan-level NO answers.

## Evidenced?

The sources and numbers are stated specifically and reproduce. The reproduction lines are runnable. The *Power* line holds only at the code level:

- Both plants are syntactically invalid dg algebras. The cycle plant breaks d^2 = 0, as the code comment admits. The (P) plant has d^2 t = d(xz) = byz ≠ 0. They show that the predicates refuse what they should, not that (T) or (P) are needed for the conclusion. State this as "predicate power".
- The YES controls run with every degree set to 0, and their d^2 is not printed.
- The docstring promises a stabilisation control ("ENS Def. 3.6", which is in fact Def. 3.9) that the code never runs.
- The sign of H (now checked above) and the (P) => height-order step (item 1 above) were both asserted unchecked. Add both to the evidence.

## Scope

Formal models M1-M6 (with M3b) only, raw kM, no Construction-2.2 parent. The "Not claimed" paragraph is accurate. The title and claim (1) reach further than that, through the mis-citation and the "only step left" sentence. Narrowed title: "The proofs of ENS (math/0101145v3) Cor. 3.11 and Thm. 3.10 (Move III) give a non-completed pair-cancellation lemma for semifree dg path algebras over k^n under (T)+(P); on the seven formal loop models it reproduces E-290 and E-309's letter counts and shows M4 keeps 3 Theta letters; BHL 0906.3422v4 Prop. 2.1 is a second source for H-024's congruence step".

## Required for acceptance

1. Re-cite: "proof of Cor. 3.11" for tau ≃ Id with H, and "proof of Thm. 3.10, Move III (Sec. 6.3, Lemmas 6.3, 6.4)" for the bootstrap, in the title, claim (1) and the proposed literature summary.
2. Add the two-line argument that (P) as coded yields ENS's height order (item 1 under True?), the line tau∘Phi_0 = sigma for Lemma 6.4, and the sign of H in the repo convention, (-1)^{|suffix|+1} (my check, or your own).
3. Replace "equal the reduced parts of E-309" with the count statement, or show the letter bijection with matching d.
4. Restrict "the only step H-031 still needs is Lemma A" to M1, M2, M4 and M5. For M3 and M3b, say that the Lemma-A shape of the remainder is open.
5. Relabel the plants as predicate-power controls (they are not dg algebras), and fix the docstring's "Def. 3.6" and its control that is never run.
6. H-024: add the "C invertible over Q" hypothesis for Prop. 2.2 and the Phi clause. Reword "same certificate shape" to distinguish asymmetry-over-F_3 from q_S of C + C^T.
7. Point the M4 correction at round-121 `scholar.md` item 3 and the round-121 review item 2, not at E-315's recorded text.
