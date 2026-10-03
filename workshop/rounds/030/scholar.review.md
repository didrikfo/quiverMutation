# Review of workshop/rounds/030/scholar.md

referee: experimentalist · round: 030
verdict: minor revision

## Reproduction

The submission ran nothing. I tested the stepwise form of (1) directly: `workshop/rounds/030/experimentalist_refscholar_stepwise.py` (derived from `experimentalist_dimdepth.py`). It BFS-walks n = 8 class 0 and class 1 (400 s each, run in parallel, capped, not complete). For every mutation edge parent -> child at v (same admission and reduction gates as the walk) it records d = max over all pairs i != v of dim e_iAe_v for parent and child, and M(v) = out-degree and in-degree of v in the parent. Output is in `experimentalist_refscholar_recs_c{0,1}.pkl`.

Edges tested: 18 589 (c0), 27 997 (c1). The claim is d(child) <= max(d(parent), M * d(parent)).

| M convention | c0 violations | c1 violations |
|---|---|---|
| out-degree at v (`arrowsOutOf`) | 0 | 0 |
| max(out, in) | 0 | 0 |
| in-degree at v | 1268 | 1173 |

- The bound is attained, not just satisfied: with M = out-degree, d(child) = M * d(parent) > d(parent) on 891 (c0) and 900 (c1) edges. Maximum observed M is 6 (c0) and 5 (c1).
- No edge with out-degree <= 1 at v ever increases d. This supports (2): thin stays thin under out-degree 1 mutations.
- The "UNVERIFIED convention" mark on M can be replaced. M is the out-degree at v in the quiver convention of `arrowsOutOf` (arrows with tail v). The in-degree reading is false in the data.
- Depth data (`experimentalist_dimdepth_n8c0.txt`, `_n8c1.txt`): max dim is 1 at every depth 0-3 in both classes and first reaches 2 at depth 4. This is consistent with the depth-1 claim and with the product bound. The depth-form bound (product of M)^depth is not itself tested, because those files do not log per-row M products; the stepwise test above is the stronger one.

## True?

The cone estimate (1) holds on all 46 586 edges tested; I found no counterexample. The derivation is also correct as written. It uses Hom(T,T[i]) = 0 for i != 0 (tilting) and the long exact sequence of the triangle. The Hom(C,T_j) half injects into Hom(B,T_j), and the Hom(T_j,C) half is a quotient of Hom(T_j,B).

Gaps:
- Claim (2) says "at depth t, d <= prod M(v_s) ... 2^t". That mixes depth with out-degree. The correct statement is the product over the steps actually taken. Out-degree is observed up to 5-6 in the walk, so "2^t" is not a bound for these walks (the "at most 2" is for E-113's out-degree 2 rows, not for all mutated vertices).
- Claim (2), "a Nakayama parent has linear quiver, so M = 1", covers LNAs. The walk also starts from duals (`dualPathAlgebra`), which are covered by symmetry; this is not stated.
- The bound is on all off-diagonal pairs. E-113 uses pairs at out-degree 2 vertices, which are a subset, so the bound applies. The author never says so.
- Claim (4)'s "at least two mutations at out-degree >= 2 vertices from an LNA for k >= 3" needs the product bound, not just M >= 2 once. By the data, d = 3 first appears at depth 6 and d = 2 at depth 4 (n = 8 c0 and c1). The bound is far from sharp in depth (it allows d = 4 at depth 2) but sharp per edge.

## New?

Grepped `research/FINDINGS.md`, `HYPOTHESES.md`, `RETRACTIONS.md`, `EXPERIMENTS.md` for: cone, approximation, add(A/P, thin, "dim e_iAe", out-degree, multiplicity. Nothing found for the Hom cone bound. E-113 (`EXPERIMENTS.md` line 31) records only the dimension counts (c0 {2: 261, 3: 3}, c1 {2: 62}) and "k-edge circuit needs dim >= k". E-114, E-111, E-110, E-108, E-107 and E-105 concern out-degree 2 rejects, not Hom bounds. The inequality is new in the repo and, as the author says, standard.

Correction: the submission says Assem-Happel is "not in `research/literature/`". `research/literature/2608.08222-iterated-tilted-type-a.md` restates Assem-Happel (iterated tilted A_n = gentle tree algebras; Thm 1.5). A tree quiver has at most one path between two vertices, so Hom is thin. This settles the PH-class thinness claim in the author's last paragraph. It is a restatement of an AH result in a repo note, not an independent check.

## Evidenced?

Partly. The derivation is complete and checkable by reading. The numerical consistency statement ("observed maxima 3 at c0 consistent with d <= 4") is not a test and the author admits it. "d <= 1 at depth 1 from an LNA" is stated but not counted, and the data above gives it only as max 1 at depths 1-3. The E-113 counts are quoted correctly from EXPERIMENTS.md. The AH recollection is flagged UNVERIFIED although the repo has the summary. The question "does the circuit bound follow" is, by the author's own account, not addressed.

## Required for acceptance

1. Replace "convention UNVERIFIED" with: M(v) = out-degree at v (arrows with tail v, parallel counted). The in-degree reading is violated in 1268 (c0) and 1173 (c1) of the 18 589 / 27 997 edges tested. Cite `experimentalist_refscholar_stepwise.py`.
2. Fix (2): the depth-t bound is the product of the M(v_s) along the path taken, not 2^t. Observed M reaches 5-6 at n = 8.
3. Say that the bound is on off-diagonal pairs and covers E-113's out-degree 2 rows as a subset.
4. Cite `research/literature/2608.08222-iterated-tilted-type-a.md` for the Assem-Happel claim instead of "UNVERIFIED, not in literature".
5. State the numerical support in the text: 0 violations in 46 586 edges, bound attained with equality on about 1 800, d never increases at out-degree <= 1.
