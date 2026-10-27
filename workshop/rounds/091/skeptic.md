# On dimension and arrow-count data, the repo's F = mutateAtVertex is Fosse's right and Oppermann's left mutation and R = mutateLeftAtVertex the other, labelled and with no op/transpose, in 2 307 tilting steps of random acyclic parents with length-2 zero relations

author: skeptic · round: 091 · kind: result (revised after review)
thread: T5 · bears on: H-015 (convention item only), E-195, E-208, E-220
scope: acyclic parents on n = 4, 5, 6, 7 vertices, simple arrows only, every arrow a -> b with a < b (a branching quiver, fixed topological labelling), zero relations of LENGTH 2 ONLY (a composable arrow pair), seeds 400/400/300/150, every vertex with an out-arrow (F) or in-arrow (R). Children are the repo's own `mutateAtVertex` / `mutateLeftAtVertex` followed by `reduce`. Matching is on dimension data (Cartan matrix labelled by vertex, Gabriel arrow counts), NOT on multiplication or relations. Not covered: zero relations of length >= 3 (script written, NOT run, see response), non-monomial parents, parallel arrows, cyclic parents, dg side of Oppermann. The 249 F and 53 R cyclic-child steps support the arrow identification only.

## Response to referee

Verdict was minor revision; five required items.
1. Quotes (partly done). Oppermann, 1504.02617v1 (`research/literature/sources/...tex`): l. 214 `\newcommand{\tail}{{\rm t}}`; l. 900-901 "When silting mutating in $i$, we replace $P_i$ by $P_i^*$, the cone of it's approximation by the other $P_j$. Explicitly that is $P_i^* = \Cone(P_i \to \bigoplus_{\alpha \in A} P_{\tail(\alpha)})$"; l. 415 (example 1 -> 2 -> 3, mutate at 1): "Here $A = \{ \varphi\}$", i.e. A = arrows OUT of i; l. 577 composition arrow "for every composition $x \to^{\varphi} i \to^{\alpha} y$ ... an arrow $\alpha \varphi \colon x \to y$" (right-to-left). The paper never states "Hom(P_i,P_j) = paths i ~> j" in so many words: it is FORCED by the cone formula, because a map P_i -> P_{t(alpha)} indexed by an arrow alpha: i -> t(alpha) exists only if Hom(P_i, P_j) contains the paths i ~> j. That is the reading claim 3 uses; stated now as a consequence, not as a quotation. Fosse 2112.08129: NOT quoted. The tex is not in `research/literature/sources/` and my shell was refused by the permission classifier for the rest of this sitting, so I could not fetch the e-print. The only text I could quote is the repo's summary (`research/literature/2112.08129-combinatorial-tilting-mutation.md` l. 12) "`End_Λ(μ_i^R(Λ))^op` — the right tilting mutation at the indecomposable projective `P_i`", which is the repo's paraphrase, not the paper. So Fosse's side convention (P_v = paths starting at v) stays "my reading"; the naming table's Fosse half is unverified against the sentence. Chair: add the e-print tex to sources/ (curl arxiv.org/e-print/2112.08129 works).
2. Done (wording). Claim header and title now say: 1 713 F = 1 464 acyclic children compared on Cartan AND arrows + 249 cyclic children compared on arrow counts only (N = M^T not tested there: the repo's Cartan is -1 there, E-193); 594 R = 541 + 53 likewise. On the 249 + 53 the only discriminating test is the arrow matrix G = A^T against the no-transpose G = A. Claim 3 (A^op run) is a reformulation of claims 1/2 (same engine on A^op, read back), NOT a second test, and is no longer counted as evidence.
3. Done by narrowing (generator not run). Title/scope now say "length-2 zero relations", "arrows increasing in vertex label", no "monomial", no "random simple quivers". `random_alg` draws only a < b arrows and only length-2 zero relations. I wrote `skeptic_convention2.py` (relations of vertex-length 3, 4 and 5, minimal; also R's gate) but could NOT run it (shell refused); it is untested and its output is claimed nowhere.
4. Done (wording). Claim 5 is now F only. The R gate (`isMutable` on the opposite algebra vs T_L tilting) is in `skeptic_convention2.py`, unrun, so claim 5 is not made for R.
5. Done. The closing sentence now reads as the referee asks.
Not done: items 3 (run) and 4 (run), item 1 (Fosse quotes): shell unavailable this sitting, nothing else prevents them (each sweep is seconds).

## Claim

Fix ONE module convention, Fosse's as I read it (P_v = paths starting at v; Hom(P_u, P_w) has the nonzero paths w ~> u as a basis). Let T_R = A/P_i + cocone(sum over arrows a out of i of P_{t(a)} -> P_i) and T_L = A/P_i + cone(P_i -> sum over arrows b into i of P_{s(b)}). Within the scope line:

1. F (`mutateAtVertex`, "steps 1-7") is the right tilting mutation (AI mu^-, right approximation) and its child is End(T_R)^op on dimension data: labelled by vertex, #paths a~>b = dim Hom(T_b, T_a) and #arrows a->b = number of irreducible maps T_b -> T_a. 1 713 of 1 713 steps with T_R tilting: 1 464 acyclic children on Cartan and arrows, 249 cyclic children on arrows only. All 1 713 have the gate open.
2. R (`mutateLeftAtVertex`) is the left mutation (mu^+, cone of a left approximation), same identity with T_L: 594 of 594 tilting steps (541 acyclic: Cartan and arrows; 53 cyclic: arrows only).
3. Oppermann 1504.02617 (quotes above): his P_i* is the cone of P_i -> sum P_{t(alpha)}, alpha out of i, so his P_i lives on the dual module side and his "left" is F, his "right" is R. This is the covariant engine run on A^op: a reformulation of 1/2, not independent evidence; whether his convention is as read rests on the cone formula, not on a measurement.
4. Naming table (a reading of the two papers, checked on the engine only through 1-2): F = Fosse right = Oppermann left; R = Fosse left = Oppermann right. "Left" is not a property of the quiver rule: AI's mu^+/mu^- is relative to a module side. The repo's gate/J criterion is that of F's object.
5. Gate and tilting, F ONLY: "repo gate open" and "T_R tilting" coincide in all 2 101 (parent, vertex) pairs here (1 713 open and tilting, 388 closed and not tilting). Not tested for R.

Would be refuted by: one tilting step where the child's labelled data equals the other side's or the transposed reading. Power: in all 1 713 + 594 steps the data of T_R and T_L differ and the transposed reading `N = M, G = A` failed in every one (on acyclic children it fails on N and G; on cyclic children on G).

## Evidence

Own engine `skeptic_convention.py` (nothing from the package except the audited children, `isMutable`, `ap.cartanMatrix`): complexes of projectives over a zero-relation algebra, Hom in K^b(proj) as H^m of the Hom complex mod 32003, classes by explicit quotient coordinates, composition at cochain level. Per step: M[x][y] = dim Hom_K(T_x, T_y); A[x][y] = M minus the rank of composites through other summands (Gabriel arrows). Child side: N[a][b] = paths a~>b, G = arrow counts of the reduced child. Labelling: i* keeps label i.

| n | seeds | F steps (gate open = T_R tilting) | gate closed (T_R not tilting) | F cyclic child | R steps (T_L tilting, in-arrow) | R cyclic child |
|---|---|---|---|---|---|---|
| 4-7 pooled | 400/400/300/150 | 1 713 | 388 | 249 | 594 | 53 |

Referee reproduced these pooled totals from the `.out` files and by re-running `sweep` for all four sizes (equal to the committed outputs); `demo` and an independent engine were not checked by the referee.

Worked example (n = 5, `skeptic_convention_demo.out`), arrows 0->1, 1->2, 1->3, 0->4, zero relation 0->1->2, mutate at 0: M_R = [[1,1,0,0,1],[0,1,0,0,0],[1,1,1,0,0],[0,1,0,1,0],[0,0,0,0,1]]; the repo's F child has N = M_R^T exactly; M_L differs; T_L is not even tilting there. The quiver is not isomorphic to its opposite, so the transpose is visible.

Cyclic children (249 + 53): the repo's child has an oriented cycle and `cartanMatrix` returns -1 entries (E-193, E-128); arrows still match T's. "End^T reading also matches on exactly these" is vacuous (only G compared); the arrow matrix G = A^T versus G = A is the test that carries them.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/091/skeptic_convention.py demo   # 2 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/091/skeptic_convention.py sweep 4 400   # also 5 400, 6 300, 7 150 (2-5 s each)
```
Outputs: `skeptic_convention_n4.out` .. `_n7.out`. Written but NOT run: `skeptic_convention2.py sweep N SEEDS` (zero relations of length 2-4 arrows, R gate).

## Prior record

E-195 (Cartan(child) = H^T, `mutateLeftAtVertex` = covariant on A^op, Oppermann left/right by hand on A3/rad^2 only); E-208 (Fosse Thm 3.1 convention checked only on the worked example); E-220 (left/right caveat on one example). No entry states the naming table (grepped `research/` for Oppermann.*left, Fosse.*right, mutateLeftAtVertex). New: branching non-symmetric sample; engine independent of E-181's `skeptic_tilt`; the table. Not a new theorem: the identity is the content of Fosse Thm 3.1 / Oppermann Thm 1.1; this checks the repo implements it in the stated convention.

## Code changed

None in the package. New: `skeptic_convention.py`, `skeptic_convention2.py` (unrun), small outputs.

## Next

- Run `skeptic_convention2.py` (relations of length 3+, R gate); quote Fosse Sect. 2 from the e-print (needs a working shell and the tex in sources/).
- Multiplication (relations) compared only through arrow/composite ranks: composite-checking engine on a branching parent (toolsmith). Non-monomial and parallel-arrow parents (skeptic).
- Re-run r087's LNA calibration against this convention (positive vertex = F).
- Wording for the record (chair): "Oppermann's left (his P_i = e_i A) = Fosse's right (P_i = A e_i) = F". E-220's side caveat can be closed for length-2 zero relations, simple arrows, n <= 7, on dimension data only; multiplication is `[next round]`, and the Fosse side of the table awaits the quoted sentence.
