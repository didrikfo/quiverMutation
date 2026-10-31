# Necessity of Z-congruence is a published one-line argument (Ladkani, math/0610685 Cor. 3.13); no independent second invariant exists for LNAs; 10 literature-certified pairs at n = 9..12 pass the Cartan tests, 3 of them the full E-269 search

author: maverick · round: 101 · kind: literature
thread: T3/T8 · bears on: H-024, E-265, E-267, E-269, F-019, idea 22
scope: n = 9..12; radical-power pairs N(n,l) only (the literature supplies nothing else); 10 pairs (list below) plus the transitive pair N(12,3)~N(12,6), not run. The E-269 exhaustive test applies only to the 3 pairs with S of signature (n-1,1), det S != 0 (n = 10, 11, 12). Bobinski-Geiss-Skowronski, Happel-Seidel and Lenzing-Meltzer-Ruan were NOT fetched (their pairs used only through Ueda's citations).
power: for the literature half n/a. For the computations: YES-only. Every pair is certified equivalent, so a NO would have voided the test; none occurred (10 of 10 char-poly and Phi-module; 3 of 3 non-vacuous q_S; 3 of 3 E-269 searches). No NO-control is added here; NO power of q_S and of E-269 is E-265 / E-269's (cross-quipu pairs). q_S is vacuous (|det S| <= 1 or det S = 0) on 7 of the 10 pairs.

## Claim

1. **H-024 (necessity step): supported by** Ladkani, *On derived equivalences of categories of sheaves over finite posets*, arXiv:math/0610685: Cor. 3.13 (derived equivalent => Euler-form matrices 1_X^{-1}, 1_Y^{-1} Z-congruent via [F]: K_0 -> K_0), Lemma 3.14 and Cor. 3.15 (M M^{-t} conjugated by the same P). Numbers are by counting `\newtheorem[section]` environments in the LaTeX (no `\label`; Prop. 3.11 = Euler form of simples = Moebius, Def. 3.12 congruent, Cor. 3.13, Lemma 3.14, Cor. 3.15). **Caveat:** the statements are for incidence algebras of posets. The proof uses only "a triangle equivalence induces an isometry of K_0 for the Euler form" and <[P_x],[S_y]> = delta_xy, so it holds for any finite-dimensional algebra of finite global dimension, LNAs included (acyclic quiver); this extension is ours, not a stated theorem. The text notes "if M1, M2 are congruent so are M1^{-1}, M2^{-1}", which is the passage from Euler matrix to Cartan matrix C. The simultaneous P with P Phi_A = Phi_B P is the proof of Lemma 3.14 (Phi is, up to sign, M M^{-t}, conjugated by the same P), not its statement.
2. **Weaker, independent statement: supported by** Bessenrodt-Holm, arXiv:math/0504547, Theorem 2.2 (first theorem of section 2, after Example 2.1 by count; label `unimod`): for any finite-dimensional algebra the unimodular *equivalence* class of C_A (invariant factors, det) is derived invariant; the paper attributes it to Bocian-Skowronski, J. reine angew. Math. 580 (2005), proof of Prop. 1.5 (UNVERIFIED, not fetched). This needs no finite global dimension but is weaker than congruence.
3. **Nikulin 1.13.3 / 1.14.2: not addressed** (no new source). E-269 does not use it; E-265's 7 q_S NOs do not use it either (q_S differing is enough for non-congruence). Only the "same genus => S isometric" upgrade would.
4. **Independent invariants (a): contradicted as candidates.** Hochschild cohomology: HH^*(A) = k for every LNA. This is F-019/idea 22's derivation, not a statement in 2312.14699; I re-checked its two inputs in the LaTeX: Bardzell resolution (label `bardzell_resn`, on-file number Thm 3.6; numbers unrecomputable, no LaTeX engine, shared `[equation]` counter) and the cochain basis p||b with p, b parallel (section 5, `\parallel{p}{b}`, s(p) = s(b), t(p) = t(b)); on a linear A_n only one path joins two vertices, and p in Gamma_n, n >= 1, lies in I, so the complex stops at degree 1. Avella-Alaminos-Geiss: **not addressed** (gentle only); F-019 shows gentle LNAs are all in the A_n class (not rechecked). Bessenrodt-Holm q-Cartan det (label `q-det`, gentle quivers, number not recomputed; formula prod (1-(-q)^k)^{c_k}): empty on LNAs, no oriented cycles gives det C_A(q) = 1 (unitriangular), the on-file conclusion, re-derived. So no independent algebraic second NO exists.
5. **What is a second NO already**: the 8 pairs of E-259 are certified inequivalent by the quipu theorem (2305.06642 via `quipuForAlmostSeparateLNA`), which is independent of Cartan data. The Cartan NOs therefore agree with an independent certificate on all 8; they are the only NOs on pairs without a quipu name that rest on H-024.
6. **(b) positive controls: supported by** Ueda, arXiv:2302.02880, Cor. 1.3 (= Cor. 5.15, label `Main 2`; intro numbering by count: Thm 1.1, 1.2, Cor 1.3, Thm 1.4, 1.5, Cor 1.6) and Cor. 1.6 (= Cor. 5.17, label `Main 4`); Ladkani, arXiv:0911.5137, Cor. 1.2 (the `c:rectline` corollary, numbered 1.2 by count and cited as [Cor 1.2] by Ueda): A(mn, m+1) ~ A(mn, n+1) ~ kA_m (x) kA_n; A(n,l) = KA_n / rad^l, digit string "l"*(n-l)+"0"*(l-2).

## Evidence

Pairs generated from the corollaries (n, l1 ~ l2); table of what each test says. "q_S" runs `disc_iso` (E-265's code); vacuous means |det S| <= 1 or det S = 0.

| pair | source | sig S | det S | char Phi equal | Z[Phi]-iso found | q_S | E-269 search |
|---|---|---|---|---|---|---|---|
| N(9,5)~N(9,3) | Ueda 1.6 | (8,0) | 0 | yes | yes | vacuous | n/a (singular) |
| N(9,7)~N(9,6) | Ueda 1.3 | (8,0) | 0 | yes | yes | vacuous | n/a |
| N(10,3)~N(10,6) | Ladkani 1.2 | (8,0) | 0 | yes | yes | vacuous | n/a |
| N(10,8)~N(10,7) | Ueda 1.3 | (9,1) | -1 | yes | yes | vacuous | YES, 5 / 22 isometries (slack 0.5 / 30) |
| N(11,6)~N(11,3) | Ueda 1.6 | (9,1) | 0 | yes | yes | vacuous | n/a |
| N(11,9)~N(11,8) | Ueda 1.3 | (10,1) | -2 | yes | yes | YES | YES, 4 / 17 |
| N(12,3)~N(12,7) | Ladkani 1.2 | (10,2) | 1 | yes | yes | vacuous | n/a (two Salem pairs) |
| N(12,5)~N(12,4) | Ueda 1.3 | (10,2) | 5 | yes | yes | YES | n/a (two Salem pairs) |
| N(12,7)~N(12,6) | Ueda 1.3 | (10,2) | 1 | yes | yes | vacuous | n/a |
| N(12,10)~N(12,9) | Ueda 1.3 | (11,1) | -3 | yes | yes | YES | YES, 4 / 15 |

"Z[Phi]-iso found" = `module_iso_search` found a unimodular P with P Phi_a = Phi_b P (a module isomorphism, not yet a congruence). E-269's search (`test`, rounds/099) returned a unimodular P with P C_a P^T = C_b for the three (9,1)/(10,1)/(11,1) pairs at slack 0.5 and again at slack 30. Two changes made in my copy, not in 099: `pq` threshold 1.5 -> 1.05 (the Salem number of these pairs is 1.1763, 1.2304, 1.2612, so the 099 code raised IndexError on them), and `M.N = n`. Hence the 099 code is not usable unchanged for Salem numbers < 1.5; the bound formula depends on lam, so a change of this kind needs the skeptic's eye. For the 7 pairs of signature (8,0), (9,1) with det 0, or (10,2) the exhaustive method does not apply (singular S, or two indefinite directions: finite-times-Z^r needs a rank-2 argument, see notebook). For definite pairs (n = 9, 10, 8 variable) the isometry group is finite and brute force would be feasible but not run.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/101/maverick_litpairs.py   # ~1 min; 10 pairs
timeout 10m .venv/bin/python -W ignore workshop/rounds/101/maverick_litiso.py     # ~1 min; 3 E-269 runs x 2 slacks
```
The venv lacked cypari and python-flint: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
Fetched (LaTeX e-prints, scratchpad): math/0610685, 2302.02880, 0911.5137, 2312.14699, math/0504547. Summary for the chair to promote: the Ladkani Cor. 3.13 / Lemma 3.14 text above is a new quotation in research/literature/math-0610685 (the existing summary says "Corollary 3.13", now matched); the Ueda, Ladkani and Bessenrodt-Holm numbers above match the on-file summaries.

## Prior record

H-024 (open); E-265, E-267, E-269; the on-file summary of math-0610685 already cites Cor. 3.13 and 3.15 and says the proof is general, but E-265/E-267 and STATE list the necessity step as unsourced, so the match was not carried over. Hochschild and AAG routes: F-019, 2312.14699 summary (idea 22), rediscovered, not new. Controls from Ueda/Ladkani pairs: the plan asked for 1 002 same-name n = 9 pairs of E-238 (experimentalist) and the 6 n = 12 within-quipu pairs (E-269, YES); my 10 are independent of that census.

## Code changed

New: `workshop/rounds/101/maverick_litpairs.py`, `maverick_litiso.py`. Nothing else touched. `module_iso_search` and `test` imported unchanged except the monkeypatched `pq`.

## Next

- Chair / skeptic: H-024 can be decided "supported, with ours: the extension from posets to finite-gldim algebras is a one-sentence proof" (statements cited above); consider an F-entry.
- Skeptic: re-derive the E-269 bound for lam in (1.1, 1.5) (the lam-range sampling of t in [lam^-1/2, lam^1/2]) before using it on other pairs.
- Maverick (later): extend the exhaustive test to two Salem pairs (rank-2 balance) so N(12,5)~N(12,4) and N(12,7)~N(12,6) are also searched; run the transitive pair N(12,3)~N(12,6).
- Scholar: fetch Happel-Seidel and Lenzing-Meltzer-Ruan for equivalent pairs with different relation lengths (all literature pairs here are radical powers, which have the same Phi-class; they give no control for the mixed-digit strings the E-259 NOs concern).
