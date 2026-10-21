# The Coxeter-key guard does not force D(i) acyclic (E-211 rerun), and inside `mutationSearchDepthFirst` step 7 is applied only at parents with no oriented cycle

author: theorist · round: 078 · kind: rediscovery
thread: T5 · bears on: H-015, E-199, E-211
scope: one explicit n = 4 monomial example (rerun, not new, at the working tree with allowCyclic=True); a reading of search.py l. 362-470 (the guarded search only; other call sites listed, not analysed); no sweep, no classification of which cyclic D(i) fail.

## Response to referee

Verdict minor revision; all four required items done.

1. Reproduction. Confirmed: at the current working tree (uncommitted `procedure.py`, `tests/test_cyclic_warning.py`) `theorist_cycD.py` raises `CyclicDescendantsError` at line 52. New wrapper `workshop/rounds/078/theorist_cycD_rerun.py` runs the 073 script unchanged with `mutateAtVertex` bound to `allowCyclic=True`. Rerun (1.7 s, this working tree): output identical to the referee's: gate True, repo child relations {al.rb - s, s.b.c, c.b}, cox x^4+x^3+x^2+x+1 (== parent: False), with rbar.b.c added (x+1)^2(x^2+1) (True), coxeterKey None None, noCycles False. The earlier Evidence figures came from the pre-edit code (round 073 tree); the guard does not change the result because the refusal is what the override bypasses.
2. Induction dropped. The argument is direct: the l. 420 test `depth > 0 and noCycles` is on the node about to be mutated, and `pathAlg.quiver` is restored from a copy of `baseQuiver` before each call at l. 437. No induction, and no LNA start is needed (a cyclic root is simply not expanded).
3. Narrowed to "the only call in search.py (l. 437)". Unguarded call sites (no `noCycles` test): `mutation.py` l. 49 and 105 (`quiverMutationAtVertices` and the second chain function), `mutation.py` l. 41 (`mutateAtVertex` via toPathAlgebra), l. 141 (dual); `pathAlgebra.py` l. 151; `procedure.py` l. 669, 671, 697 (`mutateAtVertex` direct); `lnaMoves.py` only mentions it in docstrings (l. 243, 579), no call found. These are the exposure; I did not analyse them.
4. "Next" updated (below). Title and Scope narrowed as the referee worded it.

## Claim

(1) No: `coxeterGuard` does not force D(i) acyclic. The guard compares `_coxeterKeyOrNone` of child and base; the key is None for a cyclic quiver, a None key is let through (l. 461-468), and with a cyclic start the guard is switched off (l. 362-367). It never inspects D(i).
(2) A finite-dimensional cyclic example with a cycle meeting D(i) on which the repo's `procedure.mutateAtVertex` is applied exists: E-211 (arrows a:1->2, s:2->3, b:3->4, c:4->3; relations as, sbc, cb; i = 1; true-path Cartan det 1, gl.dim <= 4). The repo child omits the relation rbar.b.c; its true Coxeter polynomial is x^4+x^3+x^2+x+1, the parent's is (x+1)^2(x^2+1). The corrected child (rbar.b.c added) has the parent's polynomial.
(3) In the guard's domain the question is vacuous. In `mutationSearchDepthFirst` mutations are tried only under `if depth > 0 and noCycles` (l. 420), where `noCycles = not simple_cycles(baseQuiver)` (l. 377) is a test on the node's whole quiver. A node with any oriented cycle (so with a cycle meeting D(i)) is recorded but never mutated, and the start is an LNA (acyclic). The test is on the mutated node itself and the quiver is restored from a copy of `baseQuiver` before the call at l. 437, so directly (no induction, no LNA start needed): inside `mutationSearchDepthFirst` (search.py l. 420), step 7 is applied only at parents with no oriented cycle; D(i) is then acyclic and E-199 applies.
Does not claim: that E-157/E-160/E-163 are affected (E-198: their parents are acyclic); anything about the other call sites above (direct callers of `mutateAtVertex`), which have no `noCycles` test and which E-211 shows are exposed; that the repo child in (2) is the only failure pattern.
Would refute (3): a mutation applied from a node whose `baseQuiver` has a cycle; i.e. any call path into `quiverMutationAtVertex`/`procedure` from the search that bypasses l. 420. I grepped: the only call in search.py is at l. 437, inside that branch. Other call sites (unguarded, not analysed): mutation.py l. 41, 49, 105, 141; pathAlgebra.py l. 151; procedure.py l. 669, 671, 697.

## Evidence

- Rerun (working tree of round 078, via `theorist_cycD_rerun.py` = the 073 script with allowCyclic=True) of `workshop/rounds/073/theorist_cycD.py`: gate admits i = 1; repo child relations are {al.rb - s, s.b.c, c.b}; no rbar.b.c; Cartan of repo child has Coxeter polynomial x^4+x^3+x^2+x+1 vs parent (x+1)^2(x^2+1); with rbar.b.c = 0 added it equals the parent's. Output printed `coxeterKey parent / child: None None` (guard lets it through) and `noCycles = False` (search would not descend).
- The failure rests on the Coxeter polynomial alone (a source tilt of a finite-gd algebra preserves it), as in E-211.
- Reading of l. 362-377, 420, 455-470 as above.

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/078/theorist_cycD_rerun.py     # about 2 s
```
The unwrapped `workshop/rounds/073/theorist_cycD.py` now raises `CyclicDescendantsError` at line 52 on this tree (the new refusal in `procedure.mutateAtVertex`); the wrapper passes `allowCyclic=True`. The example is E-211's.

## Prior record

Answered by E-211 (counterexample and the reading of the guard), E-199 (acyclic D(i) identification), E-198 (printed-path parents acyclic), E-204 (0 cyclic parents). This round adds only a rerun and the one-line induction in (3); nothing contradicts the record. No retraction involved.

## Code changed

None (new script `theorist_cycD_rerun.py` only; no library change).

## Next

- Toolsmith/chair: T5 item "does the guard force D(i) acyclic" can be closed as NO with the domain remark (3).
- The refusal in `mutateAtVertex` already exists in the working tree (uncommitted; refuses when the full subquiver D(vertex) is cyclic, `allowCyclic=True` overrides). What remains: decide whether to keep it, and the vertex-simple-tail fix for step 7.
- Theorist: make step 7 use true tails and test on random monomial cyclic algebras whether the corrected child is derived equivalent (H_1(kM), E-211 open).
