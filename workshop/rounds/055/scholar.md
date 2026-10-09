# Aihara-Iyama 2.31/2.32 and the J-criterion are verified against the LaTeX; CHZ Cor 3.6 has no "monomial" hypothesis but needs one

author: scholar · round: 055 · kind: rediscovery (citation audit; one new small counterexample)
thread: T5 / human's special request (agenda item 0) · bears on: H-015, E-066, E-122, E-126, E-128, E-136
scope: the two LaTeX sources in `research/literature/sources/`; every citation of AI Thm 2.31/2.32 and CHZ Cor 3.6 found by grep in `research/`, GLOSSARY.md, NOTES.md, `quivermutation/`, `tests/`; AI sections 1, 2 (to Thm 2.53) and 3.1-3.7 read, section 4 not; CHZ sections 1.4, 3.1-3.6, 3.12-3.21 read, sections 2 and the proof of Prop 3.18 only skimmed. Nothing re-derived beyond what is stated. No library or test changed.

## Claim

1. Aihara-Iyama Theorem 2.31 ("any mutation of a silting subcategory is again silting") and Theorem 2.32(b) (for a tilting `M` and contravariantly finite `D`, `mu^-(M;D)` is tilting iff `M` has a right `D`-approximation `g` with `Hom(g,D)` injective) are in the paper as the record cites them, as Theorems, with those numbers. Theorem 2.50(c) is the vertex-level form and states `J_i = 0` directly. The reading "silting but not tilting iff some `J_i != 0`" (E-122, E-128) is correct and needs nothing beyond 2.31 + 2.32(b); the paper does not state it as a sentence.
2. The "J = 0" premise is whole-T, not vertex-only: 2.32's iff is about the whole mutated object `D + N`, and its proof shows `Hom(N,N[<0]) = 0` follows from `Hom(N,D[<0]) = 0`; hypotheses (M tilting, D finite, generation) are automatic for `mu^-_{P_v}(A)` in `K^b(proj A)`.
3. CHZ Cor 3.6 is stated for `kQ/I` (admissible `I`), not for `kA_n/I`, and the word "monomial" does not occur in the paper. The caveat of E-066/E-122 is right in substance: Cor 3.6's path-wise form is the socle form of Prop 3.5 only when `I` is monomial (Example 3.4's identification fails otherwise). Counterexample below. For Nakayama algebras (monomial) Cor 3.6 is an iff.

What this does not claim: that the repo's rewrite equals `End(T)` (Oppermann; still Cartan-level only), or anything about section 4 of AI.

## Evidence

Source line numbers: AI = `1009.3370v3-aihara-iyama-silting-mutation.tex`, CHZ = `2509.12983v2-pavon-chz-criterion.tex`. LaTeX numbers are one counter per section (`awk` count of theorem/prop/lemma/def/remark/example/question environments; matches the summary's 2.30, 2.31, 2.32, 2.35, 2.36, 2.41, 2.43, 2.44, 3.4-3.7 and CHZ 3.5, 3.6, 3.14, 3.18-3.20 by label and content).

Statements quoted from the source:
- AI 945-946 (Thm 2.31): "Any mutation of a silting subcategory is again a silting subcategory."
- AI 982-1001 (Thm 2.32): "Let M be a tilting subcategory of T. (b) For a contravariantly finite subcategory D of M, the following are equivalent: (i) mu^-(M;D) is tilting. (ii) Any M in M has a right D-approximation g such that Hom_T(g, D) is injective." Proof printed for (a) (lines 1003-1036), "(b) dual". Definition 2.1 (353-361): silting = `Hom(M,M[>0])=0` and `thick M = T`; tilting = `Hom(M,M[!=0])=0` and `thick M = T`.
- AI 1651-1665 (Thm 2.50): "(a) T is isomorphic to right mutation mu^-_{eA}(A)... (c) T is a tilting object in K^b(proj A) iff Hom_A(eA/eA(1-e)A, (1-e)A) = 0." For `e = e_v`, `v` loopless, `eA/eA(1-e)A = S_v`, so the condition is `Hom(S_v, (+)_{j!=v} P_j) = 0`, i.e. all `J_j = 0`.
- CHZ 1382-1395 (Prop 3.5): "Let Lambda be an artin algebra... (2) the torsion pair (filt S, S^perp) induces derived equivalence iff Phi_+(S^c) in S^c" (`Phi_+ = supp soc P`).
- CHZ 1411-1431 (Cor 3.6): "Let Lambda = kQ/I be a path algebra with relations, S in Q_0... (2) (filt S, S^perp) induces derived equivalence iff every nonzero path starting outside S can be prolonged to a nonzero path ending outside S: for all p: j~>i, j in S^c, p not in I, exists q: i~>j' with j' in S^c and pq not in I."
- CHZ 1298-1330 (Ex 3.4): "the elements of soc P_i are represented by tail-maximal paths starting in i" (tail-maximal: `p not in I`, `pa in I` for every arrow `a`). This is where a monomial `I` is used silently.

Why 2.32(b) is the whole-T statement: in the proof (AI 1003-1036, mirrored for `mu^-`), `Hom(T,T[!=0])` splits into `Hom(N,D[!=0])` (iff the injectivity), `Hom(D,N[!=0])` (automatic), `Hom(D,D[!=0])` (M tilting), `Hom(N,N[!=0])` (automatic given the first, via 2.31 for positive degrees and the exact sequence for `[<0]`). So E-128's `Hom(N,N[-1])` term (which can be nonzero when some `J_i != 0`) vanishes exactly when all `J_i = 0`; E-357's "half shown" is closed, and `Hom(T,T[-1]) != 0 iff some J_i != 0` is a consequence of the proof. The generation clause is part of Theorem 2.31's proof ("the triangle shows thick mu^+ contains thick M = T", AI 949), and finiteness of `D = add(sum of other P_j)` is automatic in `K^b(proj A)`. The condition does not depend on the approximation chosen (minimal is a retract). Iteration: every step starts from a tilting `add T`, so 2.32 applies again with `M = add T`; `End(T)` is derived equivalent to `A` by Keller (AI Prop 2.3, "algebraic triangulated category").

Counterexample for Cor 3.6 (script output, `dim I = 1`): `Q`: `a:1->2, b:2->4, c:1->3, d:3->4, e:4->5`, `I = <(ab+cd)e>`, admissible and non-monomial.

| i | supp soc P_i | ends of tail-maximal paths from i |
|---|---|---|
| 1 | {4, 5} | {5} |
| 2, 3, 4, 5 | {5} | {5} |

For `S = {4}`: Prop 3.5(2) fails (`4 in Phi_+(1)`, `J_1 = Hom(S_4,P_1)` has dim 1), Cor 3.6(2) holds. By Prop 3.5 (a theorem) the torsion pair does not induce a derived equivalence, so Cor 3.6 as printed is false for this `I`. This is the E-066 mechanism (`c = [8,6,4] + [8,10,4]`, `c*(4->9) = 0`) at n = 5; the E-066 parent itself was not rerun.

Citation table (data). Status: V verified, C corrected, W wrong. "Chair" = edit to research/ the scholar may not make.

| # | citation (file:line) | status | note / proposed new wording |
|---|---|---|---|
| 1 | AI Thm 2.31 "mutation of silting is silting": EXPERIMENTS.md:303, 319, 321; syntheses/001:37, 85, 190, 195 | V | It is a Theorem (not Prop), AI 945. Needs `D` covariantly/contravariantly finite (automatic in `K^b(proj A)`). Remove "could not be checked", "cited, not re-derived" |
| 2 | AI Thm 2.32(b) iff, no monomial hypothesis: EXPERIMENTS.md:352, 355, 772, 897, 900; literature/README.md:26; tests/test_gate_without_tilting.py:5 | V | no monomial/acyclic/Hom-finite hypothesis; `M` tilting is the only one |
| 3 | "silting not tilting iff some J_i != 0": EXPERIMENTS.md:300 (E-128 title), 319, 321, 357; syntheses/001:36 | V | reading of 2.31 + 2.32(b); the `Hom(N,N[-1])` gap of :357 is closed by the proof of 2.32 (above) |
| 4 | "J = 0 implies tilting cites AI 2.32(b), not re-derived": EXPERIMENTS.md:357 | V | also AI Thm 2.50(c) states it at a vertex; Hom(N,N[<0]) automatic |
| 5 | "2.31 as Prop or Thm UNVERIFIED": EXPERIMENTS.md:245 | C | Theorem 2.31, Theorem 2.32 |
| 6 | "repo's left mutation is AI's mu^-": EXPERIMENTS.md:357 | V | AI Thm 2.50(a): the Okuyama-Rickard complex `D_v -> P_v` (P_v in degree 1) is `mu^-_{eA}(A)`, right approximation `rad`-cover |
| 7 | "minimal left approximation of P_v is sum of P_tb": EXPERIMENTS.md:364; "minimal left add(A/P_v)-approximation": :419 | C | by AI Def 2.30 a left approximation is `M -> D` (gives `mu^+`); the `P_tb -> P_v` out-arrow cover is a **right** approximation (E-355 says so). Say "right" or state the module convention |
| 8 | "J_i = H^{-1} of the mutation cone": EXPERIMENTS.md:364 | C | `J_i = Hom(N,P_i[-1])` with `N = cone(g)[-1]`; E-355 has it right. AI's `N` has `P_v` in degree 1 |
| 9 | "Hom(T,T[<0]) = 0 of AI Thm 2.32": literature/1504.02617-quivers-for-silting-mutation.md:133 | C | that is Definition 2.1(b) (tilting); Thm 2.32 says it is equivalent to injectivity of `Hom(g,D)` |
| 10 | "Theorem 2.32: an iff for 'this mutation is a derived equivalence'": literature/README.md:26 | C | iff for **tilting**; tilting implies `End(T)` derived equivalent to `A` (Prop 2.3), the converse is not claimed. Write "iff the mutated object is tilting" |
| 11 | AI Def 2.1(b), Ex 2.2(a), Prop 2.3, Thm 2.31, 2.32(b): literature/rickard-morita-theory-derived-categories.md:9, 73-74 | V | Prop 2.3 is for algebraic triangulated categories (add the word) |
| 12 | CHZ Cor 3.6 "iff", "stated for kA_n/I": literature/README.md:37 | C | stated for `kQ/I`, admissible `I`; iff proved only when `I` is monomial. "an iff for monomial `kQ/I` (all LNAs)" |
| 13 | CHZ Cor 3.6 "monomial" wording unread: EXPERIMENTS.md:774 | C | the paper has no "monomial"; Ex 3.4 (CHZ 1298) needs it |
| 14 | "Cor 3.6 path-wise would pass step 7 on the non-monomial parent, UNVERIFIED": EXPERIMENTS.md:900 | V (mechanism) | shown on a 5-vertex algebra (table above); E-066's own parent not rerun. Drop "UNVERIFIED", cite this round |
| 15 | "CHZ (Cor 3.6) not read": syntheses/001:194-196; STATE.md:24 (T5), :12 | C | both papers now read; "unread, UNVERIFIED" can go for 2.31, 2.32, 3.6 |
| 16 | summary 2509.12983: "read 2026-09-19 from the arXiv PDF" | W | untraceable (arxiv blocked in r006-053); replaced by LaTeX provenance (summary edited) |
| 17 | summary 2509.12983: add/delete-arrows remark could reshape an LNA | W | Ex 3.8 is for `kQ` without relations (summary edited) |
| 18 | summary 1009.3370: "every piecewise hereditary algebra has property (T)" as the paper's | C | Thm 3.1 is hereditary or canonical only; the extension is ours via Happel (summary edited) |
| 19 | summary 1009.3370: `mu^-_{P_i}(A)` "the cone" | C | `N = cone[-1]`, degrees 0 and 1; shift `[1]` gives the `K^{[-1,0]}` object matching CHZ (summary edited) |
| 20 | CHZ Thm 1.4, Prop 3.5, Defs/Exs 3.1-3.4, 3.9, Lemma 3.14, Props 3.18-3.19, Cor 3.20, Question 3.21, refs [14] [17] [20] | V | quoted in summary |

## Reproduction

```
timeout 2m .venv/bin/python workshop/rounds/055/scholar_chz_nonmonomial.py     # 0.5 s
grep -n "monomial" research/literature/sources/2509.12983v2-pavon-chz-criterion.tex   # no output
```

## Prior record

E-122, E-128 and E-066 already state the reductions; this round supplies the source statements, so their "cited from memory" caveats can be lifted for items 1-4, 6 above. The non-monomial caveat of Cor 3.6 is E-066's and the round-006 flag; new here: the paper has no such word, the exact place it is needed (Ex 3.4), and the 5-vertex example. New small fact, also checkable by hand: at `S = {v}`, `v` loopless, AI Thm 2.50(c) and CHZ Prop 3.5(2) are the same condition (`v` not in `supp soc P_j` for `j != v`), on any finite dimensional algebra, so the "gate vs HRS-tilt" agreement on 2052 LNAs (summary) has a three-line proof; what is left unproved is only `End(T) = repo's rewrite`.

## Code changed

None in the library. Edited: `research/literature/1009.3370-silting-mutation.md`, `research/literature/2509.12983-chz-criterion-derived-equivalences.md` (provenance, corrections, marked "r055"). Added `workshop/rounds/055/scholar_chz_nonmonomial.py`.

## Next

- Chair: apply the table's corrections (items 1-5, 7-10, 12-13, 15) to research/ and STATE.md T5; clear the STEERING special request after the referee reads the two summaries against the source.
- Skeptic: re-run the five-vertex example by hand (it is a 6-dimensional `P_1`), and decide whether item 7 is a wording slip or a convention error in E-364/E-419.
- Toolsmith: if a socle (all-paths) gate is wanted, `Phi_+` must be computed from `soc P_i` (right-multiplication kernel), not from tail-maximal paths; the repo's `J_i` already does.
- Theorist: Hom(N,N[-1]) vs J for cyclic quivers (E-128 M3) is now a corollary of 2.32's proof; the cyclic-gate caveat (simple paths only) is untouched.
