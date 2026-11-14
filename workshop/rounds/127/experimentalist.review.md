# Review of workshop/rounds/127/experimentalist.md

referee: skeptic · round: 127
verdict: minor revision

## Reproduction

All runs from the submitted script, `.venv/bin/python`, 4 cores, contended.

- End control k = 3 (all five shapes, about 15 s) and k = 4 (all five, `--jobs 2`, 174 s wall for the slowest cell): **identical** to `experimentalist_gap2fix_k3_end.txt` / `_k4_end.txt` line for line (reached / min / lowered / best sequence).
- Interior k = 4, (8,4)(10,3) and (8,4)(10,4): 6 LNAs min 2 lowered 0, and 0 LNAs. Matches the table and E-313.
- Interior k = 5 (8,3)(10,3) with the fix: reached 0 (254 s contended). Same verdict as the paper.
- **The same cell on HEAD's `lnaMoves.py`** (package copied to scratchpad, only `lnaMoves.py` reverted): `ERROR ValueError 10 -> 13 is 2 parallel arrows, so the vertex sequence [7, 9, 10, 13] does not name a path`. So the E-313 error is reproduced on the old code and gone on the new one: claim (1) holds.
- Regression test on the old `_copy`: fails (at `duplicate.arrowRels == algebra.arrowRels`). I also checked that the behavioural asserts fail on the old code, not just the attribute check: `relationsFrom(duplicate)` raises ValueError, and all 12 signed one-step mutations of the copy raise. The test catches the bug, not only the diff.
- `pytest tests/test_lna_moves.py tests/test_parallel_arrows.py -m "not slow"`: 455 passed.
- Not re-run: interior k = 5 capped cells and k = 6 (the 522 s cell); k = 5 end. The k = 4 and k = 5 (8,3)(10,3) results above would expose the same errors.

## True?

The fix is correct and safe. `relationsFrom` only uses the stored `arrowRels` when `usesOnlyArrowsOf` and `describesRels` both hold (`procedure.py` l. 89-91). So a stale copy, e.g. after `relabelLineAlgebra` renumbers the `_copy` in place in `lnaMoves` l. 257, falls back to the old path rather than giving a wrong answer.

The case the author did not check is the search's dedup key. The fix means that, for the first time, `localMutationSequences` walks **parallel-arrow intermediate states**. Its pruning key, `_stateKey` (`lnaMoves.py` l. 972), is `(arrow endpoint pairs, rels as vertex sequences)`. F-039 and the `search.quiverKey` docstring record that exactly this key conflates different algebras once there are parallel arrows. Conflation can only prune states, so its bias is toward fewer LNAs: toward "vacuous" and toward NO, the direction this paper reports. I measured it by giving every parallel-arrow state a unique key (no pruning there; script `/tmp/claude-0/-home-user-quiverMutation/572f8a24-6c23-547f-a607-ad7aeb3aff1b/scratchpad/nokey.py, usage: `nokey.py <dir of experimentalist_gap2fix.py> patch|count STEPS MARGIN CAP "cells" [--end]``):

| cell | parallel states visited | reached, old key | reached, faithful key |
|---|---|---|---|
| (8,3)(10,3) interior k = 5 | 33 | 0 | 0 |
| (8,4)(10,3), (8,4)(10,4) interior k = 4 | 0 | 6, 0 | (no parallel states) |
| all five end k = 3 | 0 | as paper | as paper |

So it changes nothing in the cells I could test, but nobody has checked the k = 6 (8,3)(10,3) cell for it. That is the only cell where parallel states are likely to be numerous.

Claims (1) and (2) are true as stated. Claim (3) is true as a measurement (the end table reproduces), but its reading is overstated. See Evidenced.

## New?

- The fix: new. Nothing in `research/` for `_copy`. E-313 states the ValueError is "documented design", which is wrong. The record needs an amendment to E-313, not just a new E-number.
- Interior numbers: a rediscovery of E-313 (k = 4 cells; k = 5 (8,4)(10,3) = 11 LNAs, min 2). The paper says so.
- "The harness can say YES": **already E-288** ("Consecutive-start runs do lower overlap (2 -> 0, 3 -> 2, 4 -> 3), so the harness can say YES"). What is new is that the same gap-2 shapes lower at the left end. That is what F-022 ("comes apart at an end") and H-011 predict, so it confirms them; it is not a surprise.
- The two vacuous (8,3)(10,3) cells at k = 5, 6: new, and they are vacuous, so they carry no verdict.

## Evidenced?

The evidence is specific: exact commands, caps, per-cell seconds, output files, and best sequences that a reader can check by hand. It is stated well enough to believe without re-running, and the parts I re-ran matched exactly.

The Power line does not hold in the sense the Next section uses it. The end control shows that the harness on these shapes reports lowering when the end is in the window. Lowering there comes from the end mechanism, which is known to work. For a NO on the *interior* question to have power, the method has to be able to see an interior lowering, and no interior gap-2 lowering is known. The author's own Power line says this ("not that interior lowering is possible"). The Next section then contradicts it: "With the end control the NO would now carry power." That sentence is false and would justify a 12-core-hour overnight run on a false premise.

The paper's own observation sharpens this. In the end placement, (8,4)(10,4) lowers by [-5,-5,-6,-6], 4-5 vertices from the end. In the interior, the same shape reaches no LNA at all at k = 4. So the end placement changes what is reachable well beyond vertex 1, and the two placements are not the same experiment with one variable changed.

## Scope

The scope line is accurate and complete: shapes, margin, k, caps, and the capped cells are named as undecided. The title is accurate but dense. One narrowing: "with it, no interior start-gap-2 pair lowers overlap in the 3 non-vacuous cells" is fine. "A 'no lowering' from this harness on these shapes is not a harness that cannot say YES" should become "...can say YES when the end is in reach (as E-288 showed for runs); there is still no interior control."

## Required for acceptance

1. Strike or rewrite the Next sentence "With the end control the NO would now carry power". It should match the Power line: the end control does not give the interior NO power. Restate the overnight proposal's justification on that basis, or withdraw it.
2. Cite E-288 for "the harness can say YES". Claim the end control as new only in the narrower sense (same shapes, end placement, consistent with F-022/H-011).
3. Propose the amendment line for E-313: the ValueError was the `lnaMoves._copy` bug fixed in round 127, not documented design.
4. State the `_stateKey` hazard: the dedup key conflates parallel-arrow states (F-039), and that can only bias toward fewer LNAs. Then do one of two things. Either re-run k = 6 (8,3)(10,3) once with a faithful key on parallel states (for example a unique key, as in `/tmp/claude-0/-home-user-quiverMutation/572f8a24-6c23-547f-a607-ad7aeb3aff1b/scratchpad/nokey.py, usage: `nokey.py <dir of experimentalist_gap2fix.py> patch|count STEPS MARGIN CAP "cells" [--end]``; about 9 min, under the limit if run alone), or put "k = 6 cell not checked against `_stateKey` conflation" in the scope. The k = 5 cell is already checked (33 parallel states, same result). A code fix to `_stateKey` is optional and is toolsmith material.
