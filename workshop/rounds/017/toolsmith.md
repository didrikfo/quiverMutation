# `reduceAgainstPivots` is now a normal form (unit test on the E-085 pair), and no monomial cord member exists for six n = 8 LNAs with a 3 early in the sequence within 6 mutation steps

author: toolsmith · round: 017 · kind: tool + negative
thread: T6 (H-017 control) and E-085 follow-up · bears on: H-017, E-085, E-087, E-084

## Claim

(1) `arrowPaths.reduceAgainstPivots` now eliminates every pivot column, not only a pivot head, so congruent combinations have equal residues. This is the library change that E-085 asked for; the E-085 pair (a = 1>2>7, b = -1>5>6>7) now has equal residues, and `theorist_step7.py` finds the kernel at target 7 that was empty. (2) At n = 8, L = 6, `MONO=1` (monomial relation sets only, no parallel arrows, arrows >= 8, rels >= 1) gives 0 members for each of the six LNAs that carry all the non-monomial members at L <= 6 (indices 4, 9, 10, 11, 12, 13 of the sorted list: 000030, 000230, 000300, 000302, 000330, 000400). Not claimed: no proof that none exists, nothing about the other 423 LNAs at L = 6 (at L = 5, indices 0-3, 5-8, 14, 15 had no cord members at all, monomial or not; 16..428 were not walked), nothing at depth 7.
Why none: not established. The two members searched in E-087 have commutativity relations, and the MONO runs find no cord member at all; a cord is a second path between two vertices, and the only way the procedure makes one is step 1/3 composites that need a two-path relation to be admissible. That is a heuristic, untested; the test would be to look at which relation of the parent produced the cord at each member.

## Evidence

Patch: 4 lines of logic. Before: the loop broke at the first non-pivot head, so a pivot column in the tail survived. After: a non-pivot head is moved to the residue and the loop continues with the tail. Pivot rows have their pivot as smallest key, so eliminating the smallest remaining pivot term never creates a smaller one; the residue has no pivot column and is canonical modulo the ideal.

New test `test_reduce_against_pivots_is_a_normal_form` in tests/test_procedure.py: an abstract case (fails on the old code, checked by stashing the patch), and a commutative square followed by one arrow (a, b congruent, residues equal, residue has no pivot column, reduction idempotent). The real-ideal half passes on both old and new code; only the abstract half detects the defect, so the E-085 pair itself is covered by rerunning `theorist_step7.py`, not by a test.

Tests (all `-m "not slow"`): tests/test_procedure.py 16 passed; tests/test_parallel_arrows.py + test_mutation_procedure.py + test_gate_without_tilting.py 26 passed, 1 xfailed (the xfail was there before). Not run: test_deeper_probing, test_fingerprint, test_invariants, test_lna_moves, test_reflections, test_relation_algebra (they import procedure/arrowPaths).

Behaviour check: the patch changes no verdict in the walk E-087 used: n = 8, L = 5, LNA 4 non-MONO gives 300 cord members with the same (arrows, rels) histogram as the round-015 output (`toolsmith_nomono_n8_L5_lna4.txt` vs `workshop/rounds/015/toolsmith_cords_n8_L5_plan.txt`). The E-084 counts themselves are not re-run (the experimentalist's task).

MONO=1 plan, n = 8, L = 6 (both directions walked, library patched), one process per LNA, 6 in parallel:

| LNA | walk | monomial cord members |
|---|---|---|
| 4 000030 | 296 s | 0 |
| 9 000230 | 197 s | 0 |
| 10 000300 | 293 s | 0 |
| 11 000302 | 196 s | 0 |
| 12 000330 | 273 s | 0 |
| 13 000400 | 315 s | 0 |

For contrast, non-MONO at L = 6 (E-087): LNA 4 has 770 and LNA 10 has 736 cord members. Sizing: one walk is 130 s alone (E-087), about 200-320 s with six running together; a search from a member at L = 6 is 250 s and 5.7e4-6.2e4 nodes. A MONO sweep over all 429 LNAs at L = 6 would be about 429 x 130 s = 15.5 h serial, but cord members at L = 5 sit in 6 of the first 16 LNAs; LNAs 16-428 are unsized (do a L = 5 `--plan` first, 10-30 s each).

## Reproduction

```
timeout 10m .venv/bin/python -m pytest -q tests/test_procedure.py -m "not slow"          # 16 passed, 5 s
timeout 10m bash -c 'MONO=1 .venv/bin/python workshop/rounds/015/toolsmith_cords.py 8 6 1 1 4 4 --plan'   # 296 s; also 9, 10, 11, 12, 13
timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 8 5 1 1 4 4 --plan                     # 31 s, 300 members
.venv/bin/python workshop/rounds/015/theorist_step7.py                                                      # 5 s, kernel nonempty, residues equal
```
Outputs: `workshop/rounds/017/toolsmith_mono_n8_L6_lna{4,9,10,11,12,13}.txt`, `toolsmith_nomono_n8_L5_lna4.txt`.

## Prior record

E-085 asked for this fix and left it as "a request to the toolsmith". E-087 recorded "no monomial cord member found at n = 6, 7" (checked small); this extends the negative to n = 8, L = 6 for six LNAs and says it was a plan over visited quivers, not a search from a candidate. E-087's controls (8 and 9 arrows, commutativity relations) stand. The n = 9 candidates are monomial (E-063, E-076), so the control still matches them only partly: this is evidence the control cannot be made monomial by choosing another member at this depth.

## Code changed

- `quivermutation/arrowPaths.py`: `reduceAgainstPivots` full reduction (docstring says why, cites E-085).
- `tests/test_procedure.py`: one new test, appended.
Callers: `procedure.py` lines ~439 (step 7 residues) and ~528 (rank/membership test); `isInIdeal` is unchanged in meaning (empty residue iff in ideal, because the old and new residues are zero together).

## Next

- experimentalist: re-run the E-084 n = 8 class 2 walk and the n = 7, 9 rejection sets under the patched rewrite; E-084's counts were made with the old one.
- theorist: is "a cord needs a sum relation" true for the mutation procedure? If it is, a monomial control at n = 8 cannot exist and H-017's n = 9 monomial candidates would need a different test (derived from a cord-free member of the class, not a cord member).
- toolsmith (next round): an L = 5 MONO `--plan` over LNAs 16-428 in shards (10-30 s each, about 2 h total) to see whether any LNA at all carries a monomial cord member; log the producing relation for each sum cord member.
