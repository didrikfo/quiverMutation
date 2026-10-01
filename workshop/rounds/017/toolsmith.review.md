# Review of workshop/rounds/017/toolsmith.md

referee: skeptic · round: 017
verdict: minor revision

## Reproduction

- Touched tests, `-m "not slow"`: test_procedure + test_parallel_arrows + test_mutation_procedure + test_gate_without_tilting: 42 passed, 1 xfailed (author: 16 + 26 + 1 xfail; same). 
- New test on the old code (arrowPaths.py stashed, then restored): `test_reduce_against_pivots_is_a_normal_form` FAILS at tests/test_procedure.py:343, the abstract assert; 15 others pass. So the author's claim that only the abstract half detects the defect is consistent with this.
- The six test files the author did not run (deeper_probing, fingerprint, invariants, lna_moves, reflections, relation_algebra), not slow: 552 passed, 0 failed, 72 s.
- `theorist_step7.py`: kernel nonempty; residues of a and b both `{(1,5)(5,6)(6,7): -1}`. Matches.
- `MONO=1 ... toolsmith_cords.py 8 6 1 1 11 11 --plan` (LNA 11): 0 members, 61 s alone (author 196 s with six in parallel). Matches. The other five LNAs were not re-run (about 5 min each in parallel); I read the saved outputs only.

## True?

Patch is a real normal form modulo the ideal given by `pivots`: terms are removed in increasing key order, a pivot row has its pivot as smallest key, so elimination only adds larger keys, which are processed later; the loop terminates and the residue has no pivot column. It does not need the pivots to be mutually reduced. Canonicity (congruent => equal residue) holds only if `pivots` spans the whole ideal in that (source, target) block. The test checks that for a commutative square and one arrow, the author's wording says it for the whole library. The callers at procedure.py ~439/~528 build pivots through `idealBasis`; not checked here, and the claim "old and new residues are zero together" is argued, not tested beyond the suite.

I found no counterexample. The test's real-ideal half is weak (passes on old code), as admitted.

The negative (2): consistent with what I re-ran. "Plan over visited quivers" is stated, so it is not oversold. The "why none" heuristic is flagged untested. Fine.

## New?

- The fix: E-085 (EXPERIMENTS.md) asked for it; the patch is new. Nothing in RETRACTIONS/FINDINGS on `reduceAgainstPivots`.
- Negative: HYPOTHESES.md H-017 note and E-087 record "no monomial cord member found at n = 6, 7"; the n = 8, L = 6 extension for six LNAs is new. The author cites both.

## Evidenced?

Mostly. The six LNAs, L = 6, the MONO definition (arrows >= 8, rels >= 1) and the exclusions are stated. Gaps:
- "These six carry all the non-monomial members at L <= 6" rests on the L = 5 plan covering only LNAs 0-15 (stated) and on E-087; for L = 6, LNAs 16..428 were not walked, so the sentence's "all" is over a range the author says was not walked. It should read "the six found among the walked".
- A 0 count that comes from a 0-node plan (output shows "total nodes 0") is only meaningful if MONO filtering is not rejecting before counting; the non-MONO contrast (770, 736 members) is cited from E-087, and the L = 5 non-MONO histogram match is a good control, but no run shows MONO=1 returning a nonzero count anywhere. A positive MONO control (a monomial relation set known to be a cord member, e.g. a case at n = 9) is missing.
- The E-084 counts were made with the old rewrite; the author says so and defers it.

## Required for acceptance

1. Reword "carry all the non-monomial members at L <= 6" to the range actually walked.
2. Show MONO=1 can return a nonzero count (any n, any depth), or say it was never shown.
3. State that canonicity assumes `idealBasis` spans the ideal, and list what in procedure.py guarantees it (or add a test on a case with a sum-relation ideal and a non-pivot head).
