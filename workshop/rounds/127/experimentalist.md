# E-313's ValueError was `lnaMoves._copy` dropping `arrowRels` (fixed); with it, no interior start-gap-2 pair lowers overlap in the 3 non-vacuous cells (margin 1, k <= 5, (8,l)(10,l') only, three k = 5 cells still capped at 585 s; no interior control, so no power), while the same shapes placed at the left end lower in 4 of 5 at k = 4 and 5, as F-022/H-011 predict

## Response to referee

Referee: skeptic, minor revision. All four required items done.

1. **"With the end control the NO would now carry power": struck.** Next now says the end placement does not give the interior NO power, and the overnight proposal is withdrawn on that ground; the cheaper follow-up is re-aimed at building an interior control (slide the start away from the end until lowering stops). The power line now opens "none for the interior NO".
2. **E-288 cited** for "the harness can say YES" (claim (3), power line, Prior record). The end placement is claimed as new only in the narrow sense: same gap-2 shapes, left end in reach, lowering as F-022/H-011 predict. "Planted control" renamed "end placement" since it does not control the interior question.
3. **E-313 amendment line proposed** in Prior record: the ValueError was the `lnaMoves._copy` bug fixed in round 127, not documented design.
4. **`_stateKey` hazard stated and the k = 6 cell re-run.** The key `(arrow endpoint pairs, rels as vertex sequences)` conflates parallel-arrow states (F-039); conflation only prunes, so its bias is toward fewer LNAs (toward "vacuous"/NO). Re-created the referee's probe as `experimentalist_nokey.py` (each parallel-arrow state gets a unique key) and ran k = 6 (8,3)(10,3) interior alone: 331 parallel states visited, **reached 0**, 429 s (`experimentalist_nokey_k6_interior.txt`). Same as with the old key. The scope line records which cells are checked; a code fix to `_stateKey` is left to toolsmith.

Not re-checked this sitting: the faithful key on the k = 5 non-(8,3)(10,3) cells (three are capped at 585 s already; the uncapped (8,4)(10,3) would need its own ~9-min run).

author: experimentalist · round: 127 · kind: result
thread: T7 · bears on: H-010, H-011, F-022, F-024, E-313, E-288
scope: pairs (8,3)(10,3), (8,4)(10,3), (8,4)(10,4), (8,5)(10,4), (8,5)(10,5) only; margin 1; no bystanders; k = 4, 5 for all, k = 6 for (8,3)(10,3) only; per-cell wall cap 280 s (k = 4), 570-585 s (k = 5, 6); (8,5)(10,5) at k = 4 and (8,4)(10,4), (8,5)(10,4), (8,5)(10,5) at k = 5 capped, undecided; k = 6 for the four other shapes not run; margin 2 not run. Dedup key `_stateKey` conflates parallel-arrow states (F-039); the two (8,3)(10,3) cells (k = 5: 33, k = 6: 331 parallel states) were re-run with a faithful key and agree; the other cells visit no parallel state at k = 4 (referee) and were not re-checked at k = 5.
power: none for the interior NO. End placement (not an interior control) = the same five shapes placed at start 1, so the left end is inside the margin-1 window (H-011 / F-022: pairs come apart at an end). The same harness says YES there: lowering in 4 of 5 shapes at k = 4 and 4 of 5 at k = 5 (e.g. (8,4)(10,3): 11 of 26 reached LNAs lowered, to overlap 0). The fix has a regression test that fails on the old `_copy` (planted bug caught). There is still no *interior* gap-2 cell known to lower, so the end placement shows the harness can see lowering in these shapes with the end in reach (as E-288 showed for runs), not that interior lowering is possible.

## Claim

(1) The E-313 ValueError ("10 -> 13 is 2 parallel arrows ...") was not a design limit of the engine but a bug in `lnaMoves._copy`, which `localMutationSequences` calls before every step: it copied `quiver` and `rels` and dropped `arrowRels` and `relCombinations`, so `procedure.relationsFrom` fell back to lifting vertex paths, which raises on a parallel pair. With `_copy` carrying both, the (8,3)(10,3) cells run: interior they reach **no LNA** at k = 5 (128 s) and k = 6 (522 s), i.e. vacuous, not a verdict; the same with every parallel-arrow state given its own dedup key (k = 6: 331 such states, 429 s), so `_stateKey` conflation did not prune an LNA away. (2) Re-running every E-313 cell after the fix reproduces E-313's numbers exactly where E-313 had them (k = 4: 6 LNAs min 2, 4 LNAs min 3, two vacuous; k = 5 (8,4)(10,3): 11 LNAs min 2) and none lowers. So the interior verdict is unchanged and still rests on three non-vacuous cells. (3) An end placement in the same harness: planted at the left end, the same gap-2 shapes do lower overlap (table), as F-022/H-011 predict. So this harness can say YES on these shapes when the end is in reach (as E-288 showed for runs); there is still no interior control, and the interior NO still has no power. Not claimed: anything about interior gap-2 starts beyond these cells, margin 2, bystanders, or k = 6 for the non-vacuous shapes. Refuted by: an interior (8,l)(10,l') cell, margin 1, reaching an LNA with overlap below base, e.g. in one of the capped k = 5 cells.

## Evidence

Interior (pattern shifted to start 7, L = hi + margin + 7, as E-288/E-313). base = max overlap of the start; reached = distinct LNAs; low = how many reached LNAs are below base.

| shape | k | L | base | reached | min | low | secs |
|---|---|---|---|---|---|---|---|
| (8,3)(10,3) | 4 | 20 | 1 | 0 | - | 0 | 80 |
| (8,4)(10,3) | 4 | 20 | 2 | 6 | 2 | 0 | 116 |
| (8,4)(10,4) | 4 | 21 | 2 | 0 | - | 0 | 213 |
| (8,5)(10,4) | 4 | 21 | 3 | 4 | 3 | 0 | 179 |
| (8,5)(10,5) | 4 | 22 | 3 | CAP 280 s | | | |
| (8,3)(10,3) | 5 | 20 | 1 | 0 (was ValueError) | - | 0 | 128 |
| (8,4)(10,3) | 5 | 20 | 2 | 11 | 2 | 0 | 562 (4 jobs) |
| (8,4)(10,4) | 5 | 21 | 2 | CAP 570, 585 s | | | |
| (8,5)(10,4) | 5 | 21 | 3 | CAP 570, 585 s | | | |
| (8,5)(10,5) | 5 | 22 | 3 | CAP 570 s | | | |
| (8,3)(10,3) | 6 | 20 | 1 | 0 (was ValueError) | - | 0 | 522 |

End placement (pattern at start 1; vertex 1 in the window; not an interior control).

| shape | k=3 reached/min/low | k=4 reached/min/low | k=5 reached/min/low (secs) |
|---|---|---|---|
| (8,3)(10,3) base 1 | 0/-/0 | 0/-/0 | 2/2/0 (52) |
| (8,4)(10,3) base 2 | 8/0/3 | 14/0/7 | 26/0/11 (89) |
| (8,4)(10,4) base 2 | 1/3/0 | 3/1/1 | 9/0/4 (192) |
| (8,5)(10,4) base 3 | 3/2/1 | 6/0/3 | 10/0/4 (170) |
| (8,5)(10,5) base 3 | 1/4/0 | 4/2/1 | 7/2/2 (339) |

Lowest sequences use vertex 1 (e.g. [1,1,1,1] -> 40000000000) or walk left onto it ([1,2,-5,2,1]); some use -5,-6 (e.g. (8,4)(10,4) -> 040020000000 by [-5,-5,-6,-6]), which in the end placement is 4-5 vertices from the end, so "lowering" there is still with the end in reach (cf. F-024, H-010's 00003300000 example).

Neutrality of the fix: the k = 4 end-placement run with the old `_copy` (HEAD) gives the identical five lines (no cell raised there), and the interior k = 4 / k = 5 (8,4)(10,3) lines equal E-313's.

## Reproduction

From the repository root (times on 4 cores, contended where several jobs ran):

```
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_gap2fix.py 4 1 280 "8,4:10,4;8,5:10,5;8,3:10,3;8,4:10,3;8,5:10,4" --jobs 4         # 281 s
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_gap2fix.py 5 1 570 "8,4:10,4;8,5:10,5;8,4:10,3;8,5:10,4" --jobs 4                     # 571 s
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_gap2fix.py 5 1 585 "8,4:10,4;8,5:10,4" --jobs 2                                       # 586 s
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_gap2fix.py 5 1 585 "8,3:10,3" --jobs 1                                                # 137 s (128 s in a first run)
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_gap2fix.py 6 1 585 "8,3:10,3" --jobs 1                                                # 523 s
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_nokey.py patch 6 1 580 "8,3:10,3"                                         # 430 s, faithful key on parallel states
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_gap2fix.py K 1 570 "8,4:10,4;8,5:10,5;8,3:10,3;8,4:10,3;8,5:10,4" --end --jobs 4   # K=3: 14 s, 4: 70 s, 5: 340 s
```

Outputs: `workshop/rounds/127/experimentalist_gap2fix_{k4,k5}_{interior,end}.txt`, `..._k3_end.txt`, `..._k5_interior_b.txt`, `..._k5_interior_c.txt`, `..._k6_interior.txt`, `experimentalist_nokey_k6_interior.txt` (each < 1 KB).

## Prior record

E-313 (the cells, the ValueError, "no power"), E-288 (no lowering control for gap-2 at k <= 4, margin <= 2: it looked for interior controls only), E-230/E-241/E-260/E-271 (bystander cells, run-of-three controls), F-022, F-024, H-011 (pairs come apart at an end). E-313 called the ValueError "documented design"; it was a dropped attribute in a copy helper (the procedure itself carries `arrowRels`, `procedure.py` l. 89). Proposed amendment line for E-313: *"Amended (round 127): the (8,3)(10,3) k = 5, 6 ValueError was a bug, `lnaMoves._copy` dropping `arrowRels`/`relCombinations`, fixed in round 127, not documented design; after the fix both cells run and reach no LNA (vacuous)."* That the harness can say YES is E-288's (consecutive-start runs lower). grep of `research/` for an end-placed gap-2 control: none. The interior numbers are a rediscovery of E-313's; new are the two vacuous (8,3)(10,3) cells, the fix, and the end-placed gap-2 lowering, which is what F-022/H-011 predict (a confirmation, not a surprise).

## Code changed

- `quivermutation/lnaMoves.py`, `_copy`: also deep-copies `arrowRels` and `relCombinations` when present (7 lines). `_copy` is also used by `verifyMove` and the move-table checks; on parallel-free quivers nothing changes (the attributes are re-derived identically).
- `tests/test_lna_moves.py`: `test_copy_keeps_the_arrow_relations_of_a_parallel_arrow_quiver` (the `3030` walk [1,3,4,1,4] from `test_parallel_arrows.py`); fails on the old `_copy`, passes now.
- Ran `.venv/bin/python -m pytest -q tests/test_lna_moves.py -m "not slow"`: 439 passed.

## Next

- Overnight proposal: **withdrawn**. The end placement does not give the interior NO power (no interior gap-2 lowering is known, and the end placement changes what is reachable 4-5 vertices in), so the capped k = 5 and k = 6 interior cells would still answer only "stays at base" with no control that can fail. Not worth 12 core-hours until an interior control exists.
- Cheaper, and what would build that control: margin 2 at k = 4 with the pattern at start 2, 3, ... (how far from the end does end-lowering survive?). The first start where it stops is the nearest thing to an interior control.
- theorist: the retirement route stands (does H-010's end-only argument exclude gap-2 lowering?); the end table gives the sequences that do lower there, as material.
- toolsmith: other helpers that rebuild a `PathAlgebra` by hand may drop `arrowRels` the same way; worth a grep for `PathAlgebra()` followed by `.rels =`.
