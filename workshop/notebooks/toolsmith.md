# Toolsmith notebook

## Round 003
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger keyed by core word; `--jobs`, `--plan`, `--summary`, `--budget-hours`). `45` at 13 pinned in tests/test_orbits_task.py.
- Not done on purpose: a key prefilter that skips walks (E-058). Ledger unsafe under two concurrent processes on one file.

## Round 006 (T3/T8)
- `rounds/006/toolsmith_orbitclass.py N`: orbit, orbit+mirror, key partitions. n = 12 48 s ... n = 16 ~32 min.
- Trap: the shell blocks `sleep N` chains; poll with `timeout .. bash -c 'until grep -q ...'`; foreground over 120 s goes to background.

## Round 010/013/015 (T6)
- `rounds/013/toolsmith_verify.py`, `toolsmith_control.py`; `rounds/015/toolsmith_cords.py N L MINREL PER FIRST LAST [--plan]` (env BOTH, MONO, HIGH, SHORT). n = 8 control: 2 members at depth 6 (E-087). Cord members only from LNAs with a 3 early (4, 9-13).
- Costs: n = 8 depth-6 walk 130 s; n = 9 candidates ~5 ms/node.

## Round 017
- Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089), test in tests/test_procedure.py. MONO=1 L = 6: 0 monomial cord members for LNAs 4, 9-13.
- Not done: LNAs 16-428 never walked for cords.

## Round 019 (T5)
- `rounds/019/toolsmith_walk.py`: scholar_walk with `--ckpt FILE`, `--max-exp N`; resume == uninterrupted. n = 8 class 2 depth 8 under the fix: 24 316 expansions, 2 rejections, 0 key-moved. Checkpoints big (38 MB): /tmp. Depth 8 ~675 s.

## Round 022 (T5)
- The 10 E-084 `M` parents are in `rounds/014/scholar_walk_n8_c2.txt` (`M (depth, rels, v, path)`); `rounds/022/toolsmith_replay.py [--old]` rebuilds them from `path` (relations match) and replays: fixed library 10/10 keep key, Cartan assertion PASS; `--old` (head-only reduction) 0/10 key, assertion FAIL 10/10.
- `procedure.mutateAtVertex(..., checkCartan=None)` / env `QM_CHECK_CARTAN=1`, `cartanDiscrepancy`, `CartanCongruenceError`; 2 tests in tests/test_procedure.py. Check costs 15 ms/step (x12 on bare mutate, +168% on gate+mutate+reduce+key); off by default. Fails on non-tilting steps by design (E-093).
- Belief: the E-085 defect is closed per step; the assertion is a Cartan-only check and cannot see a wrong rewrite with right dimensions.

## Next
- Depth 9 of n = 8 class 2 (frontier 38 907, ~35 min) as OVERNIGHT proposal; `--budget-hours` loop wrapper.
- Walk with QM_CHECK_CARTAN=1 on tilting steps only; measure whole-walk overhead (not measured).
- L = 5 MONO `--plan` over LNAs 16-428 in shards.
