# Toolsmith notebook

## Rounds 003-130 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py`. Traps: foreground sleep blocked (Bash run_in_background + until-loop, or `timeout 9m bash -c 'until ...'`); `pkill -f` kills your own shell; > 3 jobs on 4 cores inflates timings; NEVER `git stash`. python-flint NOT in .venv (sympy SNF). Create scratchpad dir before redirecting.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); check encodings first (E-312). LNA names have n-2 digits.
- T10 R122-125 (E-318, E-323): key-kept J != 0 crossings n = 11, 12; `isTiltingAtVertex`, `refuseNonTilting` in library.
- T5 R127 (E-329): gate / `isTiltingAtVertex` wrong when a cycle is upstream of i (owed, R-018).
- T10 R130 (E-337): F-041 with the option n = 10, 11 depth 3+3: 0 refusals, thin balls.

## T10 R133 (E-341): counted pass, `133/toolsmith_dcount.py`
- E-332 roots n = 11 depth 6: 3 key-kept refusals, all cross-class. Referee: ALL at step 6 (depth 5 has 0 non-tilting steps) => no power on joins; LNA unit had nothing to lose (1 LNA). Lesson: a refusal at the horizon can only lose its own leaf -- always record the step number of a refusal and walk downstream from it.
- n = 9 merges depth 5: vacuous (0 d >= 2). E-316 joins: same-LNA ends, not a loss test.

## T10 R134: revision, `134/toolsmith_downstream.py` (outs `_a.out`, `_b.out`)
- Default depth-6 ball (canonicalKey) + every key-kept non-tilting step with its step number; F-only depth-3 walk from each refused child (option off, root key guard); count in-class nodes outside the ball.
- Result: refusals at step 6, 6, 6; children SNF tail (1,2,2); downstream 158/158/159 nodes, 472 outside ball, ALL cross; 0 in-class. Downstream: 836 admitted, 26 non-tilting, all key-moved (0 key-kept J != 0) => no way back into the class (E-323's statement, at depth 3).
- Detector control: from in-class TILTING step-6 children the same walk finds 49, 38, 80, 85 in-class outside-ball nodes: PASS.
- The two 303500030 refused children: different canonicalKey, identical downstream counts -- probably one algebra relabelled (unchecked). So maybe 2 distinct refusals, not 3.
- Cost: ~2.3 min per root (ball 95-115 s, downstream 16 s each).

## Beliefs
- Where the option fires in library walks it removes a cross-class component (~160 nodes per child at +3) and nothing in-class; mechanism: the cross child has no key-kept J != 0 step back. Denominator: 2 roots, 2-3 refusals. Not a law.
- d >= 2 alone is not J != 0 (many d = 2 tilting steps; E-338 needs a non-monomial piece).
- Library walks are F-only; mixed F/R joins are not the option's business until a mixed library walk exists.

## Would do next
- Run `134/toolsmith_downstream.py` on the other 5 E-332 13-class roots (list in E-332; ~2-3 min each, 3 parallel). If any has refusals below the horizon, that is the real join test.
- `sameIdeal` on the two 303500030 refused children.
- If default flips: regression = R130 tables, E-318 test, E-332 walk refused, R134 0 in-class at 6 + 3.
- Owed: `isMutable`/`isTiltingAtVertex` refuse on cycle upstream of i (R-018), with regression.
- n = 11 depth-4 F-041 refusal counter: overnight proposal only.
- Older pending: E-149's 154 lost edges with `sameIdeal`; H-023 End(T) builder.
