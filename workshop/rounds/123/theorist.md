# (LF) holds for every graded Construction-2.2 parent (proof). On 10 small LNAs (n = 4-6, all 82 vertex instances, both sides) (LF), (H1w) and (H3) are machine-checked on real parents. Round 122's (H1) is NOT implied: its first clause fails on 26 gate-closed LNA instances and on one gate-admitted J = 0 tree step. The proof never uses that clause.

author: theorist · round: 123 · kind: proof
thread: T5 · bears on: H-023, H-031, E-319, E-235, E-256, E-252, E-279
scope: (LF): a proof for every Construction-2.2 parent that has a positive weight making d homogeneous, with all A-arrows of equal weight (any vertex i with no degree-0 loop; loops of degree >= 1 and cyclic quivers are allowed). (H1w), (H3) and H(I): machine-checked on 10 LNAs only (A4 20, 30; A5 220, 230, 320, 020, 202; A6 2220, 3300, 2400), at every vertex with an out-arrow, on A and on A^op (82 instances), plus one hand-built non-LNA tree. Parents are built by Construction 2.2 (Tate-style, exact over Q), and H = A is machine-checked on all 20. Steps are matched to the repo's F / R through E-252 (F = Oppermann's left = Def. 3.2 on A). (H3) is NOT proved for LNAs in general: g is found per instance. Hom/End(T) is not touched. No loop parent: for these no-loop steps cancellation is already E-235's result.
power: [MM] one dropped term in one parent d -> H = A fails (1 failure), NO. [LF] a real parent with positive weights but unequal A-weights (w(a2) = 5) -> NO (cycle weight -3, a letter of weight -3), which is the true answer; the round-122 toy weight-0 cycle -> NO. [SH] planted g (largest X value set to 0) on A6 3300 F1 and on the tree -> NO on both ("d0 != 0 on c2_3", "on c"). [ENG] the closed ideal (a*) in place of I -> H != 0 (2, 2 and 6 bad cells on A5 220 F1, A6 3300 F1 and the tree), NO. Strict (H1) says NO on 27 instances. That NO is true of the clause, but cancellation holds on all 27, so it is not power for cancellation.

## Claim

Setting and notation are those of E-319: kM is Def. 3.2 at i, A = the degree-0 arrows out of i, Theta = letters with a suffix, I = (Theta, d Theta).

**Proposition (LF) (proved).** Let the parent be a Construction-2.2 model of a basic finite-dimensional algebra (so d p != 0 for every letter of degree >= 1, and no word of d p is a single arrow), with no degree-0 loop at i. Give it positive letter weights w that make d homogeneous, and assume **(C)**: w(p) > w(b) for every non-A letter p with source i and every b in A. Then every oriented cycle of the kM letter quiver has positive weight, i.e. (LF) holds. (C) holds whenever all A-arrows have the same weight, in particular for the path-length grading of a homogeneous (e.g. monomial) algebra, and whenever |A| = 1. **(LF) is not implied by positive weights alone:** kQ/(a1 c) with a1: i->1, c: 1->2, a2: i->2, w = (1, 1, 5) gives f a2^-1 of weight -3 and a negative cycle f·a2*. Reweighting w(a2) = 1 repairs it.

**Correction to E-319's (H1).** Lemmas 1-2 use only **(H1w)**: every Theta letter with a linear term has a pivot x_theta (coefficient != 0, x_theta not a Theta letter, x_theta in no other Theta letter's linear part, pivots distinct). E-319's first clause, "only Theta letters have a linear term", is never used. With (H1w) in its place the Theorem (LF)+(H1)+(H3) => kM -> N q.i. keeps its proof.

**On real parents (computed).** On all 82 LNA instances (48 of them gate-admitted J = 0 steps under E-252's reading), (LF), (H1w), [PI] (pi is a dg map, d_H^2 = 0) and (H3) all hold, and the engine gives H(I) = 0. The letter quiver of kM has no cycle at all on a line. (H1w) holds at every vertex of every LNA, by a short proof (below). Strict (H1) fails on 26 instances, every one gate-closed under E-252's reading: a non-Theta letter p has d p = (a·y) + ..., where y is a letter ending at i. On LNAs, J = 0 excludes this (proof below). **It is not excluded off the line:** the tree y: 1->2, a: 2->3, a2: 2->4, z: 3->5 with relations y a and a z is gate-admitted at 2 with J = 0, and d p = a.y is linear with p non-Theta. There (H1w), [PI], (H3) (g(y) = g(p) = 1, g(c) = 2; leading word a2*·q for q/a2) and H(I) = 0 all hold. **Answer to the assignment:** (LF) is implied (given (C)). The first clause of (H1) is not implied, and it is not needed. (H3) is checked, not implied: it is the remaining open hypothesis. "Other relations through the pivot b" are harmless: A5 220 F1 (pivot x2, relation x2x3 through it, d_H c2_3 = r1·x3) and A6 2220 F1 (three overlapping relations, pivots x2, r2, c2_5) pass.

Not claimed: (H3) for every LNA; anything about loop parents of finite-dimensional algebras beyond (LF); Hom/End(T); the R/F matching beyond E-252. **Refuted by:** a parent meeting the Proposition's hypotheses with a cycle of weight <= 0; or a model with (LF)+(H1w)+(H3) and H(I) != 0.

## Evidence

**Proof of the Proposition.** Every kM letter is a.p b^-1 (prefix a when t(p) = i, suffix b when s(p) = i), with weight w(p) + w(a) - w(b), or a star b* : t(b) -> i, with weight -w(b). (i) w(p) > w_A := w(b) for b in A, for every non-A letter p from i, by induction on degree: every word of d p starts with a letter from i, which is either an A-arrow followed by a non-empty rest of positive weight, or a non-A letter from i of weight > w_A by induction. Degree-0 letters from i are all in A. This is the step that needs d p != 0 and no single-arrow word. (ii) Hence every non-star letter has positive weight: a.p and plain p trivially, and p b^-1, a.p b^-1 because w(p) > w(b). (iii) In a cycle, a star b* ends at i, and the next letter starts at i. It cannot be a star (s(b*) = t(b) != i, since there is no degree-0 loop), and it has no suffix (a suffixed letter starts at t(b) != i). So it is a.p or p with s(p) = i, of weight >= w(p) > w(b). Pairing each star with its successor gives cycle weight = sum over pairs of (w(next) - w(b)) + sum of the unpaired non-star weights > 0. **Weakest step:** "d p != 0 and no single-arrow word" is read from Construction 2.2 (each new letter kills a nonzero class, and the relations are admissible). For the 20 parents built here it holds by construction.

**(H1w) on LNAs (proof sketch).** On a line A = {a: i -> i+1}. A Theta letter is p/a, and its linear part comes only from the T-term of d(p/a) (strip a from the words a·x of d p). The base and star terms have length >= 2, because a loop at i is impossible on a line. Linear means d p contains a·x with x a single letter, which is plain, so not Theta. On a line, the cell (s, t) holds at most one Anick chain per degree, so x determines p, and the pivots are distinct and unshared. **Strict clause 1 fails => J != 0 (LNA):** the failure needs a 2-letter word y·a in some d p, with t(y) = i. Then the last relation of the chain p ends at i+1 and starts at some s <= i-1. Its part from s to i is a nonzero path killed by a, the only arrow out of i, so J_s != 0.

**Machine checks** (`theorist_lnaparent.py`, 1 s). [MM] H(parent) = A in every cell (exact), and d^2 = 0 on parent letters and on all kM letters: 0 failures on 20 parents and 82 kM's. Per instance: Theta / pivots / Theta' (a * marks a strict-(H1) failure; "(closed)" = repo gate closed or J != 0, E-252 reading).

| LNA | instances: step:Theta/pivots/Theta' |
|---|---|
| A4 20 | F1:1/1/0, F2(closed):0/0/0*, F3:0, R2(closed):0*, R3:1/1/0, R4:0 |
| A4 30 | F1:1/0/1, F2:0, F3(closed):0, R2(closed):0, R3:0, R4:1/0/1 |
| A5 220 | F1:2/2/0, F2(closed):1/1/0*, F3(closed):0*, F4:0, R2(closed):0*, R3(closed):1/1/0*, R4:2/2/0, R5:0 |
| A5 230 | F1:2/2/0, F2(closed):1/0/1*, F3:0, F4(closed):0, R2(closed):0*, R3(closed):1/1/0, R4:0, R5:2/0/2 |
| A5 320 | F1:0, F2:1/1/0, F3(closed):0*, F4:0, R2:0, R3(closed):0*, R4:1/1/0, R5:0 |
| A5 020 | F1:0, F2:1/1/0, F3(closed):0*, F4:0, R2:0, R3(closed):0*, R4:1/1/0, R5:0 |
| A5 202 | F1:1/1/0, F2(closed):0*, F3:1/1/0, F4(closed):0*, R2(closed):0*, R3:1/1/0, R4(closed):0*, R5:1/1/0 |
| A6 2220 | F1:3/3/0, F2(closed):2/2/0*, F3(closed):1/1/0*, F4(closed):0*, F5:0, R2(closed):0*, R3(closed):1/1/0*, R4(closed):2/2/0*, R5:3/3/0, R6:0 |
| A6 3300 | F1:2/1/1, F2:1/0/1, F3(closed):0, F4(closed):0*, F5:0, R2(closed):0*, R3(closed):0, R4:1/0/1, R5:2/1/1, R6:0 |
| A6 2400 | F1:2/2/0, F2(closed):1/0/1*, F3:0, F4:0, F5(closed):0, R2(closed):0*, R3(closed):1/1/0, R4:0, R5:0, R6:2/0/2 |

All 82 instances: (LF) YES (no cycle), (H1w) YES, [PI] YES, (H3) g found wherever Theta' != 0 (10 instances, 8 of them admitted), and ENG bad cells 0, sumH 0, ideal closed. The leading word is either the first two arrows after a of the stripped relation (e.g. r1/x1 in A4 30: d_H = -x2·x3 + x1*·r1, top x2·x3) or a*·p (A5 230 F2, A6 2400 R6). Both are Lemma A shapes, checked by `h031.shape`. Every gate-admitted step here has J = 0. Under the opposite F/R reading (not E-252's), 2 admitted steps (A6 3300 R4 and F2) fail strict (H1), and both still pass (H1w)+(H3).

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/123/theorist_lnaparent.py                 # all parts, controls, tree, about 1 s
.venv/bin/python -W ignore workshop/rounds/123/theorist_lnaparent.py --rows 6:3300   # one LNA
```

## Prior record

E-235 already proves cancellation for no-loop parents, conditional on the omega-reduction and Anick's lemma. E-256 checks d_1's Lemma A shape on random parents, also conditional on the omega-reduction. So **H(I) = 0 on these LNA steps is not new**. What is new: the round-122 route reaches real parents without the omega-reduction (still per instance, through (H3)); (LF) is proved, where E-319 says "checked per model, not proved"; and E-319's (H1) is corrected to (H1w). E-279 (no loops at library steps) is consistent: the letter quiver here is acyclic. A grep of `research/` for "(LF)", "locally finite", "positive weight" and "omega-reduction" finds only E-319, E-235, E-256 and H-031's line. The Proposition also covers loop parents (Step-4 letters a.p b^-1 with p: i -> i), but no finite-dimensional loop parent was built.

## Code changed

New script `workshop/rounds/123/theorist_lnaparent.py`: a Construction-2.2 builder for LNAs, kM via the 069/106 engines, the checks of `122/theorist_h031.py` (imported read-only), and a vertex-general copy of the 106 engine's cell loop. No library change, so no tests were run.

## Next

- Keeper: amend H-031 / E-319 so that (H1) reads (H1w), and record the Proposition as T5's (LF) answer. H-023's loop joint then rests on (H3) and Hom/End(T) only.
- Theorist: prove (H3) for LNA steps (the two observed leading shapes suggest g = path length with a* weighted up). Then a finite-dimensional loop parent (E-279: a cyclic quiver) through this script with a non-line builder.
- Skeptic: re-derive the Proposition's step (i) for parents with degree >= 2 letters whose d has an A-arrow-initial word of length >= 3. Run the script on n = 7 rows (`--rows`).
