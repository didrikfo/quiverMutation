# The worked examples of 2112.08129 (Ex. 3.1; Sect. 6 steps 1-2 as printed and its endpoint reachable for n = 3, 4) and 1504.02617 (Ex. 4.1, 4.2) reproduce in the library; the cyclic Ex. 3.2 does not, and Oppermann prints no loop-of-degree->=1 example


## Response to referee

1. Done. "the 16 admitted steps" is now 17 (5 + 12) in the Claim.
2. Done. Sect. 6 row narrowed (title, Claim, Scope): steps 1-2 compared with the pictures by eye; the endpoint is reachable by a 5- resp. 12-step gate/J/key-admitted sequence found by BFS (not the paper's later steps). Stated: the BFS target test was key filter + `isoTo`, and the two wrong rectangles the referee tried were rejected by the key alone. New script `toolsmith_isoNO.py` (+ `.out`) gives `isoTo` a NO on key-EQUAL non-iso targets, n = 3 and 4: the pushed-square quiver, the start A_{2n,3}, and the rectangle with the first square zero on one path (all key equal, all "iso to rect: False"); the dropped-relations and all-zero rectangles differ in key. So `isoTo` discriminates on its own.
3. Noted, not edited (research/ untouched): FINDINGS F-002 ("adds a loop on i*") is stale for the current `procedure.py`, and the `mutateAtVertex` docstring ("allowCyclic=True restores the old behaviour") is stale for r: i ~> i. Current behaviour is `KeyError (0,(2,1,0))` at procedure.py l. 221, present since at least 3eda051 (before the E-222 refusal), per the referee's check. Chair may update F-002.
4. Done. Stated in Scope: the printed finals of Ex. 3.1, 4.1, 4.2 (`want_arr`/`want_rel` and the {2->1, 1->3} set) are hand-typed from the figures, not computed; Ex. 3.1 is a re-run of `test_paper_example_acyclic`, not new evidence.

author: toolsmith · round: 102 · kind: literature
thread: T5 / T10 · bears on: H-023, H-015, F-002, E-208, E-211, E-222, E-273
scope: LaTeX in `research/literature/sources/` only (2112.08129v2 = `2112.08129v2-fosse-...tex`, 1504.02617v1 = `1504.02617v1-oppermann-...tex`; arXiv e-print not re-fetched, the repo copies are the ones E-208 / E-273 read, so versions v2 / v1 as recorded there; journal version of Oppermann not seen). Numbers are my count of `\begin{example}` per section (Fosse: Ex. 3.1 l. 280, Ex. 3.2 l. 325, Sect. 6 l. 862; Oppermann: Ex. 2.3 l. 295, Ex. 4.1 l. 403, Ex. 4.2 l. 438), not read from a compiled PDF. Computation: right mutation `quiverMutationAtVertex` + `reducePathAlgebra`, gate `mutationIsPossibleAtVertex`, J = `perI` (dim ker g_i, copied from `rounds/033`; J = 0 iff `{}`), key guard = `search._coxeterKeyOrNone`. Sect. 6 only for n = 3, 4 and with a BFS completion, not the paper's sequence (below). Printed finals for Ex. 3.1, 4.1, 4.2 are hand-typed from the figures (Ex. 3.1 only re-runs the existing test). Sect. 6: steps 1-2 compared by eye; later steps are the repo's BFS (key filter + `isoTo`; `isoTo` NO on key-equal non-iso targets shown in `toolsmith_isoNO.py`). Stale, noted not edited: F-002 ("adds a loop") and the `allowCyclic` docstring; real behaviour is KeyError, since at least 3eda051. Pictures compared as arrow/relation sets, not as dg quivers (no degrees in the library).
power: the pipeline does say NO on a planted wrong-side input: Oppermann Ex. 4.1 and 4.2 run as repo LEFT mutation (the wrong side for these labels) are refused by the gate (J = {2:1}, resp. {3:1}) and, when the rewrite is forced, by the key guard (key (1,1,1,1) -> (1,2,2,1); (1,1,1,1,1,1) -> (1,2,3,3,2,1)). So "gate True, J = 0, key kept" below is a discriminating outcome, not a vacuous one. For the cyclic case the key is None and the guard has no power (E-223).

## Claim

Per example (exact numbering in the source; "agree" = repo output equals the printed final quiver with relations, labels as printed):

| Source, example | Verdict | What the repo does |
|---|---|---|
| 2112.08129 Ex. 3.1 (7 vertices, at 3) | **supported by** | gate True, J = {}, key kept; reduced result = printed 7-vertex quiver, 4 relations, exact (already the test `test_paper_example_acyclic`) |
| 2112.08129 Ex. 3.2 (1 <=> 2, relation beta alpha: 1 ~> 1, at 1; cyclic D(1)) | **contradicted by** (library fails to reproduce the printed result) | gate True (Thm 4.4 does not apply); `mutateAtVertex` raises `CyclicDescendantsError` (E-211/E-222 refusal); with `allowCyclic=True` it does not return a wrong answer, it **crashes with `KeyError: (0, (2, 1, 0))`** (procedure.py l. 221 via 325/360); key None, so no guard; J from `perI` is meaningless on a cycle (returned `{2: -1}`). Strict xfail `test_paper_example_with_a_cycle` is still the right pin |
| 2112.08129 Sect. 6, A_{2n,3} -> kA_2 (x) kA_n, n = 3, 4 | **supported by**, steps 1-2 as printed (by eye) and reachability of the endpoint by a 5- resp. 12-step admitted sequence; **not addressed** for the paper's full sequence | steps 1 and 2 (at 1, then 3) give the printed square + zero relation 1* ~> 5 + length-3 tail (n = 3, compared by eye to the pictures); every step admitted, J = 0, key (Coxeter coefficients) kept; a gate-admitted, J = 0, key-kept sequence of 5 (n = 3: 1 3 1 2 2) resp. 12 (n = 4: 1 3 5 1 1 2 2 2 1 3 4 4) steps ends at a quiver isomorphic (quiver + relations) to the commutative 2 x n rectangle. The sequence is the repo's BFS, not the paper's (the paper's later steps are drawn for general n; I did not follow them) |
| 1504.02617 Ex. 4.1 (A_3/rad^2, at 1) | **supported by** | repo right mutation at 1: gate True, J = {}, key kept; reduced = {2->1, 1->3}, no relation = printed `phi*`, `omega` after the dg-ideal quotient. Note: Oppermann's "left" mutation matches the repo's RIGHT (Fosse) mutation in these labels; repo `left` is refused (power row) |
| 1504.02617 Ex. 4.2 (5 vertices, relations beta alpha, gamma beta, at 2, two arrows out) | **supported by** | gate True (J = {}, key kept); reduced result has arrows {3->2, 5->2, 2->4, 1->5}, relations {1->5->2, 5->2->4}: **exactly the degree-0 part of the printed final** (beta*, delta*, s, delta alpha; `d(alpha) = delta* delta alpha`, `d(s delta^-1) = s delta*`) and is isomorphic to A^op (printed claim), not to A |
| 1504.02617 Ex. 2.3 (dg model of A_3/rad^2) and the dg parts of 4.1 / 4.2 (degrees, `x` in degree 2) | **not addressed** | the library has no degrees; only H_0 of the final dg quiver is compared |
| 1504.02617: any example with a loop of degree >= 1 at i, or cyclic D(i) | **not addressed** (none printed) | Section 4 has only Ex. 4.1, 4.2 (both acyclic, loopless). The only loops of degree >= 1 are the Ginzburg loops `l_j` (degree 2, resp. d - 1) of Sect. 5 / 6, which are theorem-level (Thm. in Sect. 5, 6, Step 4 `alpha phi beta^{-1}` for a loop `phi`, l. 718 / 848), not examples, and are not representable in an ungraded library. Consistent with E-273: loops of degree >= 1 at i are inside Thm 1.1, and the library can neither confirm nor refute it |

No example is a mismatch of an acyclic step. The one disagreement is Ex. 3.2, already F-002 (strict xfail); new here is that the `allowCyclic=True` route is a `KeyError`, not a wrong quiver, and that the printed Ex. 3.2 is also an instance of the E-211 pattern (cyclic D(i)) in which the library refuses. Does not claim: that the Ex. 3.2 printed answer is right (the final quiver is said to be isomorphic to the start; I did not verify it independently, and the repo cannot), that Sect. 6 holds for all n, or anything about the J = 0 premise beyond the 17 admitted steps here (all J = 0 by `perI`, all key-kept).

A reading flagged as MINE (not in the summary `1504.02617-...md`, UNVERIFIED beyond the printed text): the summary's dictionary row "step 3's cyclic case <-> Oppermann step 4" is half right. Fosse's arrow `alpha r-bar: i* -> t(alpha)` (Ex. 3.2: `alpha beta-alpha-bar: 1* -> 2`) is, in Oppermann, **Step 2** with `phi` the degree-1 loop `omega` at i (`||alpha omega|| = ||omega|| - 1 = 0`, loops keep their degree in Step 1, l. 328); Step 4 gives in addition `alpha omega alpha^{-1}: 2 -> 2` of degree 1. A check would be to run Ex. 3.2 as the dg quiver by hand.

## Evidence

`toolsmith_litex.py` prints, per example, start, key, gate, J, raw and reduced result, key after, and the comparison with the printed final (output kept in `toolsmith_litex.out`). Totals: 3 acyclic printed finals compared (Ex. 3.1, 4.1, 4.2), 3 of 3 equal; Sect. 6: 2 steps compared with the pictures by eye (n = 3, 4 share them), 17 steps replayed with gate/J/key (5 + 12), 17 of 17 gate True, J = {}, key kept; 2 wrong-side runs refused by gate and guard; 1 cyclic example refused / crashing. Coxeter key of A_{2n,3} equals that of kA_2 (x) kA_n (n = 3, 4), as it must.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/102/toolsmith_litex.py     # 32 s, prints the table above
```

## Prior record

E-208 (Fosse: Ex. 3.1, 3.2, Sect. 6 numbered; Sect. 6 "compared by eye", Ex. 3.2 = F-002), E-199 / E-211 / E-222 (acyclic D(i) kernel = Oppermann's; cyclic D(i) refused), E-214 / E-273 (Oppermann Thm 1.1, loops of degree >= 1 allowed). The `test_mutation_procedure.py` already pins Ex. 3.1 and the Ex. 3.2 xfail. E-183 (EXPERIMENTS.md l. 685) mentions Opp. Ex. 4.2 only as "mutates at vertex 2 where step 7 is empty on both sides", its final answer "not recomputed"; I find no entry that runs Oppermann Ex. 4.1 / 4.2 through the library, nor Sect. 6 (n = 3, 4) with gate / J / key per step, nor the `KeyError` of `allowCyclic=True`. New: those three, and the Opp-left = repo-right observation (a second check of the convention that E-208 left to "the worked example").

## Code changed

None (new script `workshop/rounds/102/toolsmith_litex.py` and its output `toolsmith_litex.out`). No tests run, no library file touched.

## Next

- theorist: do Ex. 3.2 by the dg rule (Oppermann Steps 1-4 with `omega` of degree 1) by hand and compare with Fosse's printed 1*-quiver; that is the one cyclic worked example, and decides whether the printed answer is right.
- toolsmith (small, with test): make `mutateAtVertex(allowCyclic=True)` fail with a clear error instead of `KeyError` (or implement Fosse step 3 cyclic as `alpha r-bar` per arrow); only if asked, since the search never reaches it.
- chair: STATE T5 line "loops of degree >= 1 at i" can add: no printed example in 1504.02617 (v1) has one.
