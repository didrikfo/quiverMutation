# Toolsmith notebook

## Round 003
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger keyed by core word; `--jobs`, `--plan`, `--summary`, `--budget-hours`). `45` at 13 pinned in tests/test_orbits_task.py.
- Not done on purpose: a key prefilter that skips walks (E-058). Ledger unsafe under two concurrent processes on one file.

## Round 006 (T3/T8)
- `rounds/006/toolsmith_orbitclass.py N`: orbit, orbit+mirror, key partitions. n = 12 48 s ... n = 16 ~32 min.
- Trap: the shell blocks `sleep N` chains; poll with `timeout .. bash -c 'until grep -q ...'`; foreground over 120 s goes to background.

## Round 010/013/015 (T6)
- `rounds/013/toolsmith_verify.py`, `toolsmith_control.py`; `rounds/015/toolsmith_cords.py N L MINREL PER FIRST LAST [--plan]` (env BOTH, MONO, HIGH, SHORT). n = 8 control: 2 members at depth 6, not 5 (E-087). Cord members only from LNAs with a 3 early (4, 9-13).
- Costs: n = 8 depth-6 walk 130 s; n = 9 candidates ~5 ms/node.

## Round 017
- Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089), test in tests/test_procedure.py. MONO=1 L = 6: 0 monomial cord members for LNAs 4, 9-13. Heuristic (untested): a cord needs a sum relation.
- Not done: other test files importing procedure/arrowPaths; LNAs 16-428 never walked for cords.

## Round 019 (T5)
- `rounds/019/toolsmith_walk.py`: scholar_walk with `--ckpt FILE`, `--max-exp N`; exit 2 on spent slice. Resume == uninterrupted (diff, n = 7 classes 0, 1). Checkpoints are big (38 MB at n = 8 depth 8): /tmp.
- n = 8 class 2 depth 8 under the fixed library: 24 316 expansions, 63 221 algebras, 2 rejections, 0 key-moved. First 20 899 expansions vs E-084: guard-tilt 89 189 = 89 179 + the 10 key-moved. So the 10 were the E-085 defect. E-084's depth 8 was itself truncated (20 899 of 24 316).
- Slice cost: depth 7 done at 240 s, depth 8 at ~675 s total; 480 s slices fit.
- Lesson: when a budget stop printed no summary I lost a run; summary now always prints.

## Next
- Depth 9 of n = 8 class 2 (frontier 38 907, ~35 min) as an OVERNIGHT proposal; add `--budget-hours` loop wrapper.
- L = 5 MONO `--plan` over LNAs 16-428 in shards; record producing relation of each cord.
- Unit test on a real E-085 algebra.
