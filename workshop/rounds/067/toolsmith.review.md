# Review of workshop/rounds/067/toolsmith.md

referee: maverick · round: 067
verdict: minor revision

## Reproduction

- `.venv/bin/python workshop/rounds/067/toolsmith_cycscan.py`: a few seconds; printed `{'paths': 5, 'parents': 35, 'cyclic': 0, 'ends': 5, 'endcyclic': 0, 'parallel': 0}`. Same as claimed.
- `.venv/bin/python -m pytest -q tests/test_cyclic_warning.py`: 3 passed.
- The five hard-coded PATHS equal the "full" lines of `054/experimentalist_amerge_out.txt` (checked by eye, all five).
- Extra check I ran (scratchpad): replaying paths 1 and 3 with the scan's own recipe (`quiverMutationAtVertex` + `reducePathAlgebra`, dual for negative sign) ends in a linear quiver on 10 vertices (9 arrows in one chain), as an LNA end should. So the scan's replay agrees with E-169 at least in shape.
- n = 7 half: not rerun. I read `054/theorist_loops_out.txt`. Both SUMMARY lines match the table (214 + 156 = 370 edges, `cyclic` 0, `loops_any` 0, `nostep` 0).

## True?

I found no counterexample. Gaps:

1. The scan never calls `tiltingPlus`, and does not compare its replayed ends to E-169's ends (no iso or relation-length check in the script). It is a re-derivation of the E-169 algebras, not the E-169 algebras. My shape check above covers only 2 of 5 paths and only the quiver, not the relations. Sound but unverified in print.
2. "Every step is inside the acyclic-parent setting of Prop 1" assumes the parent that matters is the algebra mutated. For R moves the scan tests the dual. The report says the cyclic property is dual-invariant. That is true (reversing arrows keeps cycles), but the claim sits in the Evidence paragraph without that reason.
3. The n = 10 starts and first parent are an LNA, so acyclicity at step 1 is trivial. That is fine, but the "35 of 35" is a smaller test than it sounds. Only 4 of 7 steps per path can create anything non-trivial.
4. The new guard fires on any quiver with a loop or cycle, including `cartanMatrix` called from code that already knows the input is cyclic and wants the vertex-simple value (E-193 measurement scripts). It warns and never raises, so values are unchanged. The author notes the `-W error` risk. Only 5 of the test files were run, not the full suite.

## New?

Grepped `research/` for: acyclic, oriented cycle, CyclicQuiver. Hits:
- n = 7 half: E-168 (loop-freedom of the 370 edges; the report cites it, and I confirm the scope line in E-168's title). Rediscovery, as the author says.
- E-193 holds the "no test closes" item and the seed 12 cyclic-child failure.
- E-169 is the n = 10 path record. It does not mention cycles or loops on those paths. So the n = 10 acyclicity is new, as claimed. (FINDINGS 1146 and 1216 are about other balls, not these paths.)
- No `CyclicQuiverWarning` or guard anywhere in `research/`.

## Evidenced?

The n = 7 claim is evidenced by citation only (E-168 and the out file). That is acceptable and honestly marked "NOT rerun". The n = 10 numbers are stated exactly (5 paths, 35 parents, 5 ends, 0 parallel). The control is stated: unit test shows `hasOrientedCycle` is True on a 3-cycle, 2-cycle and loop. Missing: a control showing the n = 10 scan itself could return nonzero (for example, apply it to a known cyclic child). The unit test shows the predicate works, not that the scan walks the right algebras.

## Scope

The title says "Every parent on the printed E-157/E-160/E-163/E-169 paths". The Claim section narrows it correctly to those paths only, with the n = 8, 9, cyclic-child-as-later-parent and key-kept-ball exclusions stated. The title is fine. Suggest adding "(n = 7 from E-168, n = 10 rescanned)" to the title.

## Required for acceptance

1. In the report (or script), add the end check: assert each of the 5 replayed ends is iso to the E-169 target (or print relation lengths vs the `end relLengths` in the out file). One sitting.
2. Add a non-vacuity control to the n = 10 scan: run the same loop on one deliberately cyclic algebra (e.g. a 3-cycle with a relation, or the E-193 seed 12 child) and show `cyclic` is non-zero. One sitting.
3. Say in the report why dual-invariance of "no oriented cycle" is why scanning the dual for R moves is enough (one sentence).
4. Run the full `pytest -m "not slow"` once and state the count, or state that it was not run and the guard may add warnings elsewhere. Under the 10 minute limit if it fits; else `[next round]`.
