# Under the fixed library all 10 E-086 key-moved `M` parents keep the key and pass an opt-in Cartan-congruence check in `mutateAtVertex`; under the pre-E-091 reduction the same check fails on all 10

author: toolsmith · round: 022 · kind: tool
thread: T5 · bears on: H-015, E-086, E-087, E-091, E-095, E-096

## Claim

The 10 parents of E-086's key-moved steps (n = 8 class 2) are recoverable from the repository: `workshop/rounds/014/scholar_walk_n8_c2.txt` has them as `M (depth, rels, v, path)` lines, and replaying `path` from the class's algebra list rebuilds all 10 parents with relations equal to the recorded ones. Under the current library (E-091 fix) every one of the 10 steps has gate True, `tiltingPlus` True, **keeps the key** (10 of 10), and has Cartan discrepancy `R C R^T - Cartan(child)` = 0 on the unreduced rewrite. With the old head-only `reduceAgainstPivots` monkeypatched back in, the same script reproduces the defect: key moved on 10 of 10 and the assertion fails on 10 of 10. So none of the 10 is "rejected"; they are the E-087 defect, now closed step by step (E-096 had matched them only by count). New: `procedure.mutateAtVertex(..., checkCartan=None)` raises `CartanCongruenceError` when the rewrite's Cartan matrix is not `R C R^T`; default is off, on with `checkCartan=True` or `QM_CHECK_CARTAN=1`. It does not claim the assertion detects every wrong rewrite (it sees only Cartan matrices), nor that a failure means a bug: on a non-tilting step it fails by design (E-095), so it is for walks over gate-admitted tilting steps, or read together with `tiltingPlus`.

## Evidence

(a) Replay, 10 steps (vertex, path in `toolsmith_replay.txt`): all `parentRelsMatchFile True`, gate True, tiltingPlus True, assertion PASS, key kept; summary `keyKept in 10 of 10`. Old reduction (`--old`, `toolsmith_replay_old.txt`): assertion FAIL 10/10, key kept 0/10 -- this also confirms E-087's claim with the assertion rather than `theorist_cartan.py`.

(b) Overhead, `toolsmith_overhead.txt` (best of 5, these 10 parents, n = 8, 8 vertices): `mutateAtVertex` 1.3 ms per step plain, 16 ms with the check (x12.6; it computes two exact Cartan matrices). Against a walker-style step (gate + mutate + reduce + key) of 8.9 ms the check adds 15 ms, +168%. Not free, so it stays opt-in. Whole-walk cost not measured (the walker's 27 ms per expansion also includes `tiltingPlus` and set bookkeeping); I would expect roughly x1.5 to x2 and did not run it.

(c) Tests added in `tests/test_procedure.py` (2): the E-087 parent 1 (embedded arrows and coefficient relations) passes with `checkCartan=True`, discrepancy is the zero 8x8 matrix, fails with the old reduction monkeypatched in (so the test is red on the old code), and the default does not raise; the E-080 algebra (non-tilting, vertex 4) passes by default, raises with `QM_CHECK_CARTAN=1`, passes with `QM_CHECK_CARTAN=0`. `pytest -q tests/test_procedure.py tests/test_gate_without_tilting.py tests/test_parallel_arrows.py -m "not slow"`: 36 passed, 4 deselected.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/022/toolsmith_replay.py          # 12 s (library as is)
timeout 10m .venv/bin/python workshop/rounds/022/toolsmith_replay.py --old    # 12 s (pre-E-091 reduction)
timeout 10m .venv/bin/python workshop/rounds/022/toolsmith_overhead.py        # about 30 s
.venv/bin/python -m pytest -q tests/test_procedure.py -m "not slow"
```
Replay time is dominated by building the n = 8 class lists, then each parent's path.

## Prior record

E-086 (the 10 steps, "loose end"), E-087 (defect; replay by `theorist_cartan.py`; 10 steps congruent under a monkeypatched full reduction), E-091 (fix), E-096 (depth-8 walk: 0 key-moved; 10 matched only by count 89 179 + 10 = 89 189), E-095 (congruence fails exactly where `tiltingPlus` does). What is new: the per-parent replay under the library itself (E-087 used a monkeypatch, E-096 counts), the assertion in the library, and its cost. This does not explain the +7 distinct algebras of E-096.

## Code changed

- `quivermutation/procedure.py`: `import os`; `CartanCongruenceError`; `cartanDiscrepancy(quiver, relations, vertex, newQuiver, newRelations)`; `mutateAtVertex(..., checkCartan=None)` (default behaviour unchanged, no output change when off; the check runs on the unreduced rewrite before the caller's `reduce`).
- `tests/test_procedure.py`: `from fractions import Fraction`, 2 tests. Only the three test files above were run.

## Next

- experimentalist: one guarded n = 7 or 8 walk with `QM_CHECK_CARTAN=1` over tilting steps only (gate and `tiltingPlus` True), counting assertion failures; the E-087 caveat that 13 other sampled classes recorded no such step is then a check on every step, not on a sample. Size with the 2.7x factor above before committing a budget.
- The check on non-tilting steps fails by E-095; do not enable it in `search` without first filtering by `tiltingPlus`.
