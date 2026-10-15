# The uniform transpose in E-181 is the `cartanMatrix` indexing (verified for LNAs and duals n <= 7 at dimension and arrow level, plus `mutateLeftAtVertex` = the opposite-algebra construction); identification with Oppermann's left/right is by hand on A3/rad^2 only

author: skeptic · round: 066 · kind: result (explains the E-181 uniform transpose; does NOT close the left/right item; mostly a rediscovery, one new quiver-level check)
thread: T5 · bears on: H-015, E-181, E-189, E-193
scope: hand rebuild of ONE example (A3/rad^2 = Oppermann Example 2.3/APR); sweep of all LNAs and duals n = 4..7 (30+112+420+1584 steps, every vertex with an out-arrow, no gate, J = 0 and J != 0 alike), acyclic, monomial. Same H code as E-181 (a re-run, n <= 7 of E-181's n <= 8). `mutateLeftAtVertex` checked against the opposite-algebra construction on A3/rad^2 and LNAs/duals n = 4..7 (2146 steps), not against Oppermann's printed right mutation. Not covered: parallel-arrow or non-monomial parents, cyclic parents, n >= 8, a hand A5 instance, any claim about Oppermann's proof (only his statement, his example and his notation).

## Response to referee

1. `\tail` (item 1): done. `mutation.tex` l. 214 only has `\newcommand{\tail}{{\rm t}}`; "target" is my inference from l. 321 (A = degree-0 arrows starting in i), l. 901 (P_{\tail(alpha)} in the cone) and l. 909 (`\start(phi) != i = \tail(phi)`). Claim 3 and Prior record reworded.
2. `mutateLeftAtVertex` (item 2): checked, `skeptic_left.py 4 5 6 7`, output `skeptic_left.out`, under a minute. Dual construction = the round-050 complex T built on A^op (`dualPathAlgebra`), H_B = dim Hom there. Result: Cartan(mutateLeftAtVertex(A,v)) = H_B (not H_B^T) at 2146 of 2146 LNA/dual steps n = 4..7, plus A3/rad^2 at v = 2, 3 (printed); H_B never symmetric; every arrow of the reduced left child has H_B[j][i] > 0 (0 violations in 12 098 arrows; the same-index reading fails at all). Honest limit: `mutateLeftAtVertex` is by its code `_opposite(mutateAtVertex(_opposite(.)))`, so this checks the code against an independent H on A^op, not that A^op-covariant is Oppermann's RIGHT mutation of A. That last identification (Hom(-,A) of his P_i^* cone) is still a hand argument, unchecked; also 2146-step run includes J != 0 steps (Hom(T,T[-1]) != 0 in 880 of them, as in E-181), where Cartan still matches H_B.
3. Docstring item (item 3): Next now reads "if the Oppermann-right identification holds"; not stated as settled.
4. Same H code (item 4): yes. `skeptic_side.py` execs the same `skeptic_tilt.py` complex code as `theorist_enddim.py` (lines 11-12 identical). It is a re-run of E-181 at n <= 7, so "independent" is dropped. It is independent of `procedure.py`/`tiltingPlus` only.
5. Title/"closes"/l. 404 ff. (item 5): title and kind narrowed. The hand rebuild was compared with Oppermann's APR example at `mutation.tex` l. 404 ff. only through its arrow/relation shape (2->1, 1->3, 2->3 = composite); I did not re-read l. 404 ff. line by line this sitting and do not claim to match his printed text. Also dropped the arrow-level count as a separate result (nearly forced once C = H^T).

## Claim

Three statements, each of which a counterexample would refute.
1. `cartanMatrix` stores C[j][i] = number of paths i~>j mod I (code: `matrix[target][source]`; docstring "entry (j,i) is dim e_j A e_i", functional composition). With P_x the covariant projective, Hom(P_x,P_y) = paths x~>y, so C = H^T *by definition at T = A*. The "uniform H^T" of E-181 is therefore the same indexing applied to the child, `Cartan(child) = H(T)^T`, not a Cartan misreading. It would only be a misreading if the child's C were H, which is false at 100% of steps below.
2. The side: the repo's End(T) is End of the covariant projectives, T_k = P_k -> (+_{k->h} P_h) (maps along arrows OUT of k). In functional composition (e_jAe_i = paths i->j, Oppermann's notation) these are right modules e_xA; in the repo's traversal-order tuples (`p + (b,)` = p then b) the same objects are left modules A e_x. Same category, two names; do not read "left modules" in `skeptic_tilt.py` and "right dg modules" as a conflict. The quiver of the child is End(T), not End(T)^op: every arrow i->j of the reduced child has H[i][j] > 0, and the transposed reading (H[j][i] > 0) fails at every arrow.
3. Oppermann's Def. 3.2 reading in E-189/E-193 is on the same side. His A = arrows of degree 0 out of i, P_i^* = Cone(P_i -> (+_{alpha in A} P_{t(alpha)})) (`mutation.tex` l. 901). `\tail` is only a macro for `t` (l. 214); "target" is inferred from l. 321, 901, 909. His APR example (l. 404 ff.) has the same arrow/relation shape as the repo's output (below; compared by shape, not line by line). The theorist's q-prefix model (E-193: q a path starting at i, red drops q nonempty) is the traversal-order image of his "rightmost arrow in A" (functional order).

Does NOT claim: that the repo's rewrite presents End(T) (G2 beyond the acyclic span statement, the cancellation step of Thm 1.1 and J = 0 are untouched); that nothing could be mirrored in a non-LNA case.

## Evidence

Hand rebuild (A3/rad^2: phi: 1->2, psi: 2->3, psi.phi = 0; Oppermann writes composition right to left, the repo's relation is the tuple ((1,2),(2,3))). T = T_1 + P_2 + P_3, T_1 = (P_1 --phi--> P_2) in degrees (-1,0). By hand: Hom(P_2,T_1) = 1 (identity on the degree-0 term; a chain map T_1 -> P_2 needs lambda.phi = 0, so lambda = 0, hence Hom(T_1,P_2) = 0); Hom(T_1,P_3) = <psi> (psi.phi = 0 makes it a chain map); Hom(T_2,T_3) = <psi>; T_3 and T_1 not related the other way; no Hom in degrees +-1. So End(T) has the arrow 2->1 (phi*) and 1->3 (omega), 2->3 equals the composite: the tilted A3 Oppermann draws. Script output (`skeptic_side.out`): H = [[1,0,1],[1,1,1],[0,0,1]], Hom(T,T[+-1]) = 0, repo child (unreduced) arrows (1,3),(2,1),(2,3) with relation (2,1,3) = (2,3) (his d(omega phi^-1) = omega phi* - psi), Cartan(child) = [[1,1,0],[0,1,0],[1,1,1]] = H^T and != H. Parent: Cartan = [[1,0,0],[1,1,0],[0,1,1]] = paths i~>j in row j.

Sweep (sum over n, all steps, the round-050 complex code, the same H as E-181):

| n | steps | H not symmetric | C = H^T | C = H only | neither | child arrows | arrows with H[i][j] = 0 (repo reading) | arrows with H[j][i] = 0 (transposed reading) |
|---|---|---|---|---|---|---|---|---|
| 4 | 30 | 30 | 30 | 0 | 0 | 84 | 0 | 84 |
| 5 | 112 | 112 | 112 | 0 | 0 | 434 | 0 | 434 |
| 6 | 420 | 420 | 420 | 0 | 0 | 2076 | 0 | 2076 |
| 7 | 1584 | 1584 | 1584 | 0 | 0 | 9504 | 0 | 9504 |

Power: H is never symmetric, so C = H^T and C = H are never both true; the test could fail at every step. Arrows are those of `procedure.reduce(child)`: in the unreduced child 8 (n = 4) arrows lie under a length-1 relation and are zero, so they are not End(T) arrows (first version of the script counted them; fixed). The arrow column is a quiver-level statement E-181 did not make, but it is nearly forced once C = H^T holds (acyclic End(T)), so it is a second view of the same fact, not a separate result.

Theorist's open item (E-193 request): "if every step outputs End(T)^op in the convention where the parent is End(P)^op each step is still a derived equivalence; a mix would not be". The data exclude a mix at dimension and quiver level for every LNA step n <= 7; the opposite-everywhere reading is the same set of objects renamed: the expected duality (Hom(-,A) takes T to a right-mutation object with End(T)^op; `mutateLeftAtVertex` = opposite, mutate, opposite) is an ARGUMENT, not a result here: `skeptic_left.py` shows only that `mutateLeftAtVertex` agrees with the covariant construction on A^op (see Response, item 2). Whether the library's "left" is Oppermann's right is a conjecture.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/066/skeptic_side.py 4 5 6 7     # 12 s, output in skeptic_side.out
```

## Prior record

E-181 already states the uniform transpose ("one uniform op convention, never a mix") and E-189/E-193 list "the left/right convention" as unclosed. `research/` has no entry that says why the transpose is uniform. New here: (a) it is the `cartanMatrix` indexing, an identity at T = A, so nothing in E-181 needed a convention to be chosen; (b) the quiver-level arrow check; (c) the usage-based reading of `\tail` (l. 321/901/909) and the hand rebuild of the A3 APR example, matching by shape; (d) `mutateLeftAtVertex` = the A^op construction at Cartan level, n <= 7. `RETRACTIONS.md` grep for "transpose" finds nothing relevant.

## Code changed

None (new script and output only). No tests run, no library file touched.

## Next

- Toolsmith: the trap is documentation: add one line to `cartanMatrix`/`procedure.mutateAtVertex` docstrings (draft wording, conditional). Not done by me. Only if the Oppermann-right identification is checked (item 2 limit); until then say only "`mutateAtVertex` acts on covariant projectives; `mutateLeftAtVertex` is the same on the opposite algebra".
- Theorist/scholar: one hand A5 instance of Def. 3.2 against `procedure.py` 340-366 (still undone; this note only checks the A3 example and the side).
- The side is explained for LNAs n <= 7 and the one hand example; the left/right NAMING is not settled; non-monomial or parallel parents were not rerun in the quiver-level form (E-181 covers 38 parallel steps at dimension level).
