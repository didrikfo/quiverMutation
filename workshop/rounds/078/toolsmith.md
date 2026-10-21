# All 8 E-216 chains (4 children, 4 parents; n = 7 classes 1, 2) have recorded F/R move lists, are Hom-verified edge by edge, and child and parent end at the same seed in 3 of the 4 reps

author: toolsmith · round: 078 · kind: result + tool
thread: T10 · bears on: E-213, E-216, E-212, E-211
scope: n = 7; one representative per E-213 presentation class (class 1 once, class 2 three times), child and parent = 8 algebras, not the 14 edges; target ball depth 4 (complete), algebra ball depth 5 (stops at first level with a hit); generation of K^b(proj) by T assumed (E-168); child and parent end at the same seed in reps 0, 1, 3 but at different seeds (1 vs 8) in rep 2, whose seed-to-seed link is still open.

## Response to referee

1. Broken recorded scripts: fixed minimally. `workshop/rounds/073/theorist_cycD.py` (line 52) and `workshop/rounds/065/theorist_g2cyclic.py` (lines 89, 119) now pass `allowCyclic=True` (comment `r078`). Outputs checked against the pre-guard code (HEAD worktree, scripts as recorded): byte-identical for both (12 and 6 lines), so E-211's reproduction command runs again with unchanged output. Grep of `workshop/` for other direct callers (`mutateAtVertex`, `mutateLeftAtVertex`, `mutateAtVertices`): 12 other scripts; 7 ran to the end with no `CyclicDescendantsError` (066/skeptic_left, 066/skeptic_side, 065/theorist_g2exact, 065/theorist_g2cartan, 015/theorist_step7, 022/toolsmith_replay, 022/toolsmith_overhead); 5 need arguments and were not run (034/skeptic_outdeg3, 031/toolsmith_parallel, 059/theorist_enddim, 059/theorist_enddim2, 027/scholar_pairtest), so not checked.
2. Tests: added for `mutateLeftAtVertex`, `mutateAtVertices` (positive and negative vertices) and loops (at the vertex, below it; loop above accepted). 6 pass in `tests/test_cyclic_warning.py`; full non-slow suite 1429 passed, 1 xfailed.
3. Seed index: printed per chain (`seed index i of N`); outputs `toolsmith_chains_r0..3_out.txt` regenerated, move lists and verdicts identical to before, 8 of 8 VERIFIED. The seed ends are in the Evidence section; the first draft's "seed-to-seed link not addressed" is replaced by the actual result, which includes a rep (2) whose child and parent end at different seeds.
4. Retitled "Hom-verified"; "verified" means Hom-consistent (Hom(T,T[+-1]) = 0 and Hom dims = Cartan) plus sign-aware structure iso at the meet, not a proof of derived equivalence. "Total edges 32" now stated as a sum.

## Claim

The J = 0 iso-meet chains of E-216 are now explicit move lists, and each of the 8 verifies: every edge has Hom(T,T[-1]) = Hom(T,T[1]) = 0 and dim Hom(T_i,T_j) = Cartan of the child (E-212 `horn`, independent of `tiltingPlus`, gate, key), and each meet is a sign-aware isomorphism (and passes the weaker E-212 structure-graph `iso`). 8 of 8 verify, 0 fail. Seed indices (new, printed in the output): chains end at seeds of `classes[base]` as follows (rep 0: index 1 of 12; rep 1: 8 of 14 for both; rep 2: child 1, parent 8 of 14; rep 3: 1 of 14 for both). Reps 0, 1, 3 share the seed between child and parent; rep 2 does NOT (child ends at seed 1, parent at seed 8), so for rep 2 the link between seeds 1 and 8 is still open. Not claimed: that child ~ parent for rep 2, equivalence beyond the Hom test, that the search is minimal (depths are the first hit, as E-216), or generation. A refuted claim would be any edge below printing `!` (Hom nonzero) or `NOSTEP`; none occurs. The reps are the E-216 ones (parsed from `074/skeptic_failkids_an_out.txt`); the chains found here are a re-search, with the same totals as E-216 (2, 2, 2, 3 / 6, 6, 5, 6).

## Evidence

Moves `kind v`: F = tilt at v of the algebra, R = same on the opposite algebra. Alg side is applied to the child/parent; seed side is applied to a seed of the class (a seed is an LNA or dual) and meets the alg side up to isomorphism.

| rep (class, v) | side | alg-side moves | seed-side moves | edges | Hom + Cartan | meet |
|---|---|---|---|---|---|---|
| 0 (1, 7) | child | (none: child ~ node) | F7 F5 | 2 | pass | iso |
| 0 | parent | F4 F7 | F7 F5 F7 F5 | 6 | pass | iso |
| 1 (2, 3) | child | - | R6 R5 | 2 | pass | iso |
| 1 | parent | F4 F3 | F1 F1 F2 F4 | 6 | pass | iso |
| 2 (2, 6) | child | - | F6 F7 | 2 | pass | iso |
| 2 | parent | F5 | F1 F1 F2 R7 | 5 | pass | iso |
| 3 (2, 7) | child | - | F7 F6 R3 | 3 | pass | iso |
| 3 | parent | F3 F7 | F6 F7 F5 F5 | 6 | pass | iso |

Seed ends (index in `classes[base]`): rep 0 child 1 / parent 1; rep 1 child 8 / parent 8; rep 2 child 1 / parent 8; rep 3 child 1 / parent 1. Edge count is the sum of the per-chain counts 2+6+2+6+2+5+3+6 = 32 (a sum of the printed `edges`, not an independent count), all Hom-verified. Raw output: `toolsmith_chains_r0_out.txt` .. `r3_out.txt` (per chain: moves, per-edge verdicts, meet verdicts).

Power caveat (E-212): the verifier's NO is thin (1 Hom rejection in 195 corruptions); a pass here is consistent with, not proof of, derived equivalence of each chain's ends beyond the Hom test plus E-168.

### Code change (E-211 refusal)
`quivermutation/procedure.py`: `mutateAtVertex(..., allowCyclic=False)` now raises `CyclicDescendantsError` (a `ValueError`) when D(vertex) (full subquiver on vertex and its descendants) has an oriented cycle or loop. A cycle elsewhere is not refused. `allowCyclic=True` restores the old behaviour. Test: `tests/test_cyclic_warning.py::test_mutateAtVertex_refuses_cyclic_descendants` (the E-211 example refused at 1, accepted with the flag; cycle above the vertex accepted). Tests run: `test_cyclic_warning, test_procedure, test_mutation_procedure, test_parallel_arrows, test_cycles, test_search_soundness, test_gate_without_tilting, test_double_mutation, test_edge_moves` with `-m "not slow"`: all pass. Behaviour change: direct callers (`mutation.quiverMutationAtVertex`, `mutateAtVertices`, `mutateLeftAtVertex`) that mutate at a vertex with cyclic D(i) now raise; the search never does (E-204: 0 cyclic parents, acyclic only). Slow tests were not run; the full non-slow suite was run after the revision: 1429 passed, 1 xfailed. Chair to decide whether to keep. Extra tests added after review: `test_left_and_multi_mutation_refuse_cyclic_descendants` (mutateAtVertices with negative vertex = left, mutateLeftAtVertex, the guard on the opposite quiver; left at 1 with nothing above is accepted) and `test_loop_at_or_below_the_vertex_refused` (loop below the vertex, loop at the vertex, loop above only accepted).

## Reproduction

```
mkdir -p /tmp/claude-0/ts75
timeout 10m .venv/bin/python -W ignore workshop/rounds/078/toolsmith_chains.py R 4 5 > workshop/rounds/078/toolsmith_chains_rR_out.txt   # R = 0..3, ~60-85 s target ball + <5 s; four together ~110 s
timeout 10m .venv/bin/python -m pytest -q tests/test_cyclic_warning.py -m "not slow"
timeout 10m .venv/bin/python -W ignore workshop/rounds/073/theorist_cycD.py     # now with allowCyclic=True; same output as before the guard
timeout 10m .venv/bin/python -W ignore workshop/rounds/065/theorist_g2cyclic.py
```

## Prior record

E-216 (depths only, "move lists not recorded"), E-213, E-212 (Hom verifier, class 2 records 0-8), E-211. No RETRACTIONS entry bears on this. New content: recorded chains and their Hom verification at the 4 unjoined classes; the E-211 guard.
