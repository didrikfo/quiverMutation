# Every parent on the printed E-157/E-160/E-163/E-169 paths is acyclic (n = 7 from E-168, not rerun; n = 10 rescanned, ends checked against E-169), and cyclic input now warns

author: toolsmith · round: 067 · kind: rediscovery + tool
thread: T10/T5 · bears on: E-193, E-168, E-169, H-015
scope: n = 7 classes 1, 2 (44 paths, 370 edges, from E-168, NOT rerun) and n = 10 group-A merge (5 paths, 35 edges, rerun; ends iso to the E-169 targets 5/5, non-vacuity control 4/5 on injected cycle); the algebras along the paths only, not the whole key-kept ball

## Response to referee

Verdict taken: minor revision. All four required items done; scripts `workshop/rounds/067/toolsmith_cycscan.py` (unchanged) and new `toolsmith_cycscan2.py` (`.venv/bin/python workshop/rounds/067/toolsmith_cycscan2.py`, seconds).

1. End check (done). Each of the 5 replayed ends is compared with the LNA built from the `end relLengths` printed in `054/experimentalist_amerge_out.txt` (`[4,5,0,5,5,0,0,0]`, `[5,5,5,0,4,4,0,0]`, `[6,0,5,0,4,0,3,0]` x3), as a structure graph (vertices, arrows, relation paths, via `procedure.relationsFrom`). Result: 5 of 5 isomorphic (relation counts 4, 5, 4, 4, 4 on both sides); 0 of 5 labelled-equal, which is expected, since the replayed vertex labels differ from the LNA's. So the scan replays the E-169 algebras up to isomorphism at the end, on all five paths, with relations (not only quiver shape, not only 2 of 5). The intermediate parents are not compared to anything; they are the replay itself.
2. Non-vacuity (done). Same loop, with an extra arrow 10 -> 1 added to a copy of the step-4 parent of each path (never mutated): `cyclic` = 4 of 5 paths (35 parents scanned), so the scan can return nonzero on the algebras it walks. The fifth parent has no directed path 1 -> 10, so the added arrow closes no cycle there; I did not tune it. A bare 3-cycle quiver gives `hasOrientedCycle` True. I did not use the E-193 seed 12 child (needs that run's pickles; not this sitting).
3. Dual-invariance (done, one sentence). For an R move the scan tests the dual of the algebra, which is enough because dualising reverses every arrow, and reversing all arrows maps oriented cycles (and loops) to oriented cycles, so the algebra has one iff its dual does.
4. Full suite (done). `.venv/bin/python -m pytest -q -m "not slow"`: 1426 passed, 1928 deselected, 1 xfailed, 4 warnings, 414 s. One of the warnings is a `CyclicQuiverWarning` (from `invariants.py:109` calling `cartanMatrix`, in a test that feeds a cyclic quiver); the others are a networkx hash-change UserWarning and similar, not from the guard. No failures from the guard. The `-W error` caveat stands.

Changed in the text: title (referee's suggestion plus the end check), Evidence/Claim as below. Gap 3 (35 of 35 is a smaller test than it sounds: step 1 is an LNA) is accepted as stated: the substantive tests are the later parents; the count 35 is what was scanned, not 35 independent trials. Gap 4 (guard fires on callers that know the input is cyclic) unchanged, warns only.

Updated claim: n = 10 half: 35 of 35 parents and 5 of 5 ends without loop or oriented cycle, on replays whose ends are isomorphic to the E-169 targets (5 of 5), and the scan returns 4 of 5 on a cyclic injection. n = 7 half unchanged (cited, NOT rerun). Scope line below adds the end check.

## Claim

On the printed paths none of the 405 algebras mutated (370 at n = 7, 35 at n = 10) and none of the 49 path ends has a loop or an oriented cycle, so every step on them is inside the acyclic-parent setting of Prop 1 / E-193 (Lemma L). The n = 7 half was already recorded (E-168, round 054), so E-193's "what no test closes: every parent on a printed path is acyclic" is closed for these paths, and the chair's agenda item was stale on that point. The n = 10 half is new (35 of 35, 0 cyclic). It does NOT cover: the 25 failing steps' children (cyclic children are possible, E-193 seed 12) as parents of later steps off the printed paths; the key-kept balls; n = 8, 9. A single cyclic parent on a printed path would refute it.

A cheap guard now exists: `cartanMatrix`, `arrowPaths.isInIdeal` and `_minimalGenerators` emit a `CyclicQuiverWarning` (never raise, value unchanged) on a quiver with a loop or an oriented cycle, since `allPathsBetween` is vertex-simple and the answer can be silently wrong there (E-193).

## Evidence

| set | paths | edges (parents) | cyclic parents | cyclic ends | source |
|---|---|---|---|---|---|
| n = 7 class 1 (E-157 child+parent, E-160, E-163) | 26 | 214 | 0 | 0 | E-168, `054/theorist_loops_out.txt` |
| n = 7 class 2 | 18 | 156 | 0 | 0 | same |
| n = 10 E-169 (group-A merge) | 5 | 35 | 0 | 0 | `toolsmith_cycscan.py`, new |

The n = 7 check in `054/theorist_loops.py` tests the algebra mutated at each edge (dual for R moves; the oriented-cycle property is dual-invariant) with `nx.simple_cycles` and self-loops, and the end of each side; 0 `NOSTEP`, so it replayed every printed edge. The n = 10 scan rebuilds the five signed paths of `054/experimentalist_amerge_out.txt` (sign = dual) with no `tiltingPlus` call; it reports 0 parallel-arrow parents too. Control for non-vacuity: the injected-cycle run of `toolsmith_cycscan2.py` gives cyclic = 4 of 5 paths (see Response); the new unit test shows `hasOrientedCycle` is True for a 3-cycle, 2-cycle and a loop, and `test_cycles.py::test_the_search_ends_the_branch_at_a_cycle...` now emits the warning (a cyclic child does arise in the search).
Unit test: `tests/test_cyclic_warning.py` (acyclic input silent under `simplefilter("error")`; the three entry points warn on a 3-cycle with answer unchanged).

## Reproduction

```
.venv/bin/python workshop/rounds/067/toolsmith_cycscan.py        # n = 10, a few seconds; last line {'paths': 5, 'parents': 35, 'cyclic': 0, ...}
.venv/bin/python -m pytest -q tests/test_cyclic_warning.py       # 3 passed
cat workshop/rounds/054/theorist_loops_out.txt                   # n = 7 summary (519 s + 392 s to rebuild pickles; not rerun)
```

## Prior record

E-168 (round 054) already ran the n = 7 loop/oriented-cycle scan: "0 loops, 0 oriented cycles" on 370 edges; E-193's open list and STATE T5 repeated "every parent on a printed path is acyclic" as unclosed, so this is a rediscovery for n = 7 and E-168 should be cited there. New: n = 10 (E-169 paths, which E-168 did not cover), and the code guard. E-193 records the cyclic-Cartan failure (seed 12, i = 3); the warning is its cheap consequence. Nothing in RETRACTIONS bears on this.

## Code changed

- `quivermutation/arrowPaths.py`: `CyclicQuiverWarning`, `hasOrientedCycle`, `warnIfCyclic`; called from `cartanMatrix` and `isInIdeal` (`import warnings` added). Behaviour of values unchanged.
- `quivermutation/procedure.py`: `_minimalGenerators` calls `ap.warnIfCyclic`.
- `tests/test_cyclic_warning.py` (new). Ran: test_cyclic_warning, test_cycles, test_parallel_arrows, test_procedure, test_mutation_procedure (`-m "not slow"`): 53 passed, 1 xfailed (pre-existing), 1 warning (the intended one in test_cycles). Full `-m "not slow"` suite later run: 1426 passed, 1 xfailed, 414 s (see Response). Not guarded: `relationAlgebra.isInIdeal`/`cartanMatrixExact` (vertex model) are not guarded. The warning prints once per call site under default filters, so a search through cyclic children is not flooded, but a run with `-W error` would now fail on cyclic quivers (none seen in the tests run).

## Next

- chair/scholar: record in E-193's open list and STATE T5 that the printed-path acyclicity is closed by E-168 + this round, leaving the cancellation step of Thm 1.1 and the left/right convention.
- theorist: the open case is a cyclic *child* becoming a parent off the printed paths (key-kept balls); I can count cyclic parents in the E-192 balls (n = 7 first) if wanted.
- toolsmith (later): make the cyclic Cartan right (non-simple paths up to the nilpotency bound) only if a thread needs cyclic input.
