# Toolsmith notebook

## Rounds 003-130 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py`. Traps: foreground sleep blocked (Bash run_in_background + until-loop, or `timeout 9m bash -c 'until ...'`); `pkill -f` kills your own shell; > 3 jobs on 4 cores inflates timings; NEVER `git stash`. python-flint NOT in .venv (sympy SNF). Create scratchpad dir before redirecting.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); check encodings first (E-312). LNA names have n-2 digits.
- T10 R122-125 (E-318, E-323): key-kept J != 0 crossings n = 11, 12; `isTiltingAtVertex`, `refuseNonTilting` in library.
- T5 R127 (E-329): gate / `isTiltingAtVertex` wrong when a cycle is upstream of i (owed, R-018).
- T10 R130 (E-337): F-041 with the option n = 10, 11 depth 3+3: 0 refusals, thin balls.
- T10 R133-134 (E-341, E-344): 3 key-kept refusals at E-332's 2 crossing roots, all step 6 (= horizon), 2 algebras; F-only +3/+4 from them: 0 in-class outside ball. Lesson: record the step of every refusal; a horizon refusal loses only its leaf.

## T10 R137: decision data, `137/toolsmith_t10walk.py` (+ `.out`)
- Pre-registered rule in the submission (HARMLESS = L_in 0 and 0 cross nodes on; non-vacuous = key-kept refusal at step <= D-1; R1-R4 regressions).
- Depth 7, library F-only + Visited, 9 walks (7 roots root side + 303500030 dual + 303305000 root): L_in 0, G 0 everywhere; L_x 19 / 12 at the 2 crossing roots; on balls 0 cross nodes. Key-kept refusals: step 6: 3 (E-341's), step 7: 20 new, all cross. NO non-tilting step at steps 1-5 in any walk.
- R1 F-041 n = 8 10/10, n = 9 52/52 on = off, 0 refusals. R2 E-318, R3 E-323, R4 seeds 514 (v5 F), 3277 (v3 R): PASS.
- Planted power: refuse every step into 5 in-class step-4 nodes (depth 5) -> L_in = 9. PASS.
- Verdict by rule: "supports flipping default on, scoped" -- but non-vacuity rests on the same 3 step-6 refusals (1 step below). Expect the referee to push on that; my bar for non-vacuous was low (should have asked for >= 2 steps below).
- Cost: depth-7 walk 57-472 s (3 parallel), nodes 868-6950; depth 8 ~3x, so not in 10 min. `shard`/`merge` mode (by step-1 child, own Visited) verified = walk at D = 5 (681 keys), overlap ~1.9x.
- Timing caveat: my off walk calls `isTiltingAtVertex` at every admitted step (spy), so off/on times are equal by construction; the option's real cost is unmeasured.

## Beliefs
- In the 13-class of n = 11 group 5, the option removes only cross-class subtrees; non-tilting key-kept steps appear only at steps >= 6 from the LNA roots, and every one has a cross-class child (E-323's "no way back" holds again). Denominator: 9 walks, 23 refusals. Not a law.
- d >= 2 alone is not J != 0 (E-338 needs a non-monomial piece).
- Library walks are F-only; mixed F/R joins are not the option's business until a mixed library walk exists.

## Would do next
- If the conference wants >= 2 steps below a refusal: depth 8 via `shard ROOT SIDE 8 off|on 6 m` at the 2 crossing roots (~24 commands <= 8 min) -> OVERNIGHT proposal; depth 9 ~3x more.
- Measure the option's cost honestly: off walk WITHOUT the spy vs on walk, same root, alone on the machine.
- The 5 unrun duals (003033030, 023033030, 203033030, 223033030, 305500030 dual), ~1-8 min each.
- If default flips: regression = R130 tables, E-318 test, E-332 walk refused, R137 R1-R4; update test `test_refuse_non_tilting_drops_exactly_e318s_step_and_is_off_by_default`.
- Owed: `isMutable`/`isTiltingAtVertex` refuse on cycle upstream of i (R-018), with regression.
- Older pending: E-149's 154 lost edges with `sameIdeal`; H-023 End(T) builder.
