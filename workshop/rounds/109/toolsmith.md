# On every shortest E-287 tilting path B_n -> A_n (n = 5, 6, 7; shifts <= 7, first-depth hits, depth 3, 4, 5) the repo's `mutationIsPossibleAtVertex` admits at every step and the reduced result equals End(T') (labelled presentation; Cartan matrix and canonicalKey), although the intermediate End(T) are not LNAs

author: toolsmith · round: 109 · kind: result
thread: T10 · bears on: H-027, E-187, E-280, E-282, E-287
scope: n = 5, 6, 7 only; all SHORTEST tilting-exchange paths (4, 5, 6 paths) to the first-depth A_n hits (2 hits each), BFS over all single-summand exchanges of rigid K0-unimodular complexes among the n(n-1) exceptional objects, shifts <= 7, depth <= 5 (n = 5, 6) / 5 (n = 7, stopped at the first hit). Longer paths (depth +1 hits: 5, 7 more at n = 5, 6) not audited. Generation (Bruning) still UNVERIFIED. Not a statement about the LNA-only move set. Says nothing about the Coxeter-key guard: the key is preserved by derived equivalence on tilting steps, so it cannot refuse here.
power: the gate and the match test can say NO, shown on the control: at B_6 the gate refuses 4 of 12 (vertex, side) cells (v1L, v4L, v5R, v6R) and the repo result equals a tilting neighbour's End in exactly the 8 admitted cells (= E-282's 8; B_5: 6 of 10; B_7: 9 of 12). The match test also says NO for the wrong side: on all 62 audited steps exactly one of {right, left} reproduces End(T'). The key column is not evidence (implied by tilting); no claim about the guard.
backlog: n/a

## Response to referee

1. Title/claim reworded: the key-kept column is implied by derived equivalence, so it is no longer claimed as evidence; "guard keeps" is dropped from title, claim, scope and power. The key-refusal control (a J != 0 step; round-102 cyclic refusal) is moved to *Next*. Done.
2. Match definition stated in the claim: labelled Cartan matrix, transpose accepted in the test (`C == H2 or C == H2.T`). Upgraded for all 62 steps with new `toolsmith_relcheck.py` (outputs `toolsmith_relcheck_{5,6,7}.out`): the matching-side reduced result has the same `fingerprint.canonicalKey` (labelled arrows plus relation set) as `reducePathAlgebra` of End(T'), 12/12, 20/20, 30/30, none by the dual, none by transpose. Caveat: my End(T') presentation can list a redundant relation (n = 5: 1345 follows from 132 and the commutative square), which gave 22 spurious mismatches until End(T') was reduced by the repo's reduction; so "equal" means equal after the repo's reduction, and the referee's "equality of presentations" is met in that sense. Done.
3. Gate-True cells that fail the match: at the B_n control, none (every gate-True slot matches a tilting exchange; the 4/4/3 failing cells (n = 5/6/7) are all gate-False). Along the paths, yes, on the non-matching side: of 124 (step, side) cells, 62 are the matching side (all gate True, all match), and 38 more are gate True with a non-matching result (35 left, 3 right; 24 are gate False). These are cells the gate admits whose result is not End(T'): the gate is not a test of "is the tilting neighbour". Done.

## Claim

Let T_0 = B_n, T_1, ..., T_d (d = 3, 4, 5 for n = 5, 6, 7) be any shortest chain of tilting exchanges ending at End(T_d) = A_n. At every one of the 12 (n = 5), 20 (n = 6), 30 (n = 7) steps, the repo's single mutation at the exchanged vertex (right at 9/16/25 steps, left at 3/4/5) passes `mutationIsPossibleAtVertex`, its reduced result has exactly the labelled Cartan matrix of End(T_{k+1}) (transpose accepted in the test; the identity reading was the one that matched at all 62 steps) and, after reducing End(T')'s presentation, the same canonicalKey (same arrows and relations), (coxeterKey unchanged, but that is implied by derived equivalence and is not evidence about the guard). So no gate refusal occurs along these paths: the gate does not obstruct the multi-step tilting route, and E-280's "A_n not in orbit(B_n)" is not caused by a refusal at a step of these paths. What does obstruct is the scope of the gate: E-187 is a theorem for LNAs only, and the intermediate End(T) are not LNAs, so "gate True" there is the repo's procedure run outside the proved range. It does not claim that every tilting exchange is admitted: of the 12 vertex-sides at B_n the gate refuses 4 (n = 6), each of whose slots still has a tilting exchange (reached by the other side, or non-vertex-local). Refuted by a shortest path whose step is refused.

## Evidence

Presentation of End(T): thin (all Hom <= 1-dimensional, checked), triangular; arrows = non-factorable nonzero Hom (composites tested in homotopy classes); relations = minimal zero paths plus equality of parallel nonzero paths (sign check: no sign issue flagged). Convention fixed by End(B_n) = repo B_n (lone relation 2345 / 3456 / 4567). Per T_k: LNA? gentle?

| n | paths | steps | T_0 | T_1..T_{d-1} (count of End(T)) | T_d | gate True | key kept | repo result = End(T') |
|---|---|---|---|---|---|---|---|---|
| 5 | 4 | 12 | LNA, not gentle | 8 not LNA (5 arrows with commutativity relation: 5; 4 arrows, 2 in-arrows at one vertex: 3) | A_5 (LNA, rad^4 = 0) | 12/12 | 12/12 | 12/12 |
| 6 | 5 | 20 | same | 15 not LNA (6 arrows 9; 5 arrows 6) | A_6 | 20/20 | 20/20 | 20/20 |
| 7 | 6 | 30 | same | 24 not LNA (7 arrows 14; 6 arrows 10) | A_7 | 30/30 | 30/30 | 30/30 |

(Counts of End(T) are over path-steps, so repeated nodes count repeatedly.) Examples (n = 5, labels 1..5, repo reading): T_0 = B_5; T_1 = arrows 1->3, 2->5, 3->2, 3->4, 4->5, relations 132-path, 1345, 3-2-5 = 3-4-5 (commutative square) -- a non-gentle, non-monomial, non-line algebra; T_2 has relation 3125 and 345; T_3 = A_5 (relation of length 4, `43125`). No intermediate is a gentle LNA, and no intermediate is an LNA at all; the first step is an LNA-to-non-LNA move (a vertex-local E-282 step, End not A_n) that the gate admits. Where the gate or key "first refuses": nowhere on these paths. n = 7 depth: no hit at depth <= 4 (nodes 12, 76, 327, 1046 new at depths 1..4; 2667 at depth 5 with 2 hits). Depth-to-first-hit 3, 4, 5 at n = 5, 6, 7: one more per n, 3 data points (I do not call it a law).
Side facts: the matching side is always unique; the other side's gate is often True too but its result is a different algebra. Tilting exchange vs. repo mutation agree at every step (no step where the repo's operation is a silting-only mutation).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/107/toolsmith_cones.py 7   # 2m38s, writes the objects file; I moved it to workshop/rounds/109/toolsmith_objs_7.json
timeout 10m .venv/bin/python workshop/rounds/109/toolsmith_pathaudit.py 5 5 > workshop/rounds/109/toolsmith_pathaudit_5.out   # 3 s
timeout 10m .venv/bin/python workshop/rounds/109/toolsmith_pathaudit.py 6 5 > workshop/rounds/109/toolsmith_pathaudit_6.out   # 6 s
timeout 10m .venv/bin/python workshop/rounds/109/toolsmith_pathaudit.py 7 6 > workshop/rounds/109/toolsmith_pathaudit_7.out   # 35 s
```
timeout 10m .venv/bin/python workshop/rounds/109/toolsmith_relcheck.py 7 6   # (also 5 5, 6 5) canonicalKey comparison, < 1 min
(The cones step must be re-run, or the json moved, for n = 7; n = 5, 6 read round 107's files. The script prints the control first, then every path.)

## Prior record

E-287 (hits, paths not audited, "intermediate End not checked to be LNAs"): this answers that and adds n = 7. E-282 (one-step: gate agrees with the independent test on 281/281): consistent; my B_n control reproduces its 8 (n = 6). E-187 (gate rule proved for LNAs): the new fact is that the same procedure also admits and matches on non-LNA End. E-280: "A_n not in orbit(B_n)" computed with LNA-only moves; the present audit shows the tilting route is not blocked by gate or key, so the orbit statement reflects the move set (vertex mutations that stay in LNAs), as E-287 guessed. `RETRACTIONS.md` mentions E-280 only as the source of the correction of E-218 ("A_n, B_n not derived equivalent"); this result is consistent with that correction, not affected.

## Code changed

None in the repository. New: `toolsmith_pathaudit.py` (execs 107's tiltsearch/cones preludes; uses `tests/helpers.path_algebra`, `mutation`, `reduction`, `search._coxeterKeyOrNone`), `.out` files (< 15 KB), `toolsmith_objs_7.json`. No tests touched.

## Next

- Skeptic: Cartan-and-canonicalKey equality is done (response 2); a check up to isomorphism rather than labelled equality not needed. Also thin presentation sign check only for +-1.
- Theorist: why does the repo procedure admit non-LNA algebras with a commutative square (is E-187's step 1-7 reading valid there)? And is the depth n-2 first-hit pattern provable (A_n and B_n both in the D_n derived class)?
- Toolsmith/experimentalist: key-refusal control (moved from the claim; the guard is untested by this run: find a tilting-valid step whose mutation changes the key, or a non-tilting step) for this pipeline (find a tilting-valid step whose mutation changes the key, e.g. a J != 0 step of E-268), and n = 8, 9 (72 objects; cones would exceed 10 min, so an OVERNIGHT proposal: cones 8 then pathaudit 8 7).
