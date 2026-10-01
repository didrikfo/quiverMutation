# Toolsmith notebook

## Round 003
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger keyed by core word; `--jobs`, `--plan`, `--summary`, `--budget-hours`). `45` at 13 pinned in tests/test_orbits_task.py.
- Left out on purpose: a key prefilter that skips walks (E-058). Ledger is not safe under two concurrent processes on one file.

## Round 006 (T3/T8)
- `rounds/006/toolsmith_orbitclass.py N`: orbit, orbit+mirror, Coxeter-key partitions. Cost n = 12 48 s ... n = 16 ~32 min.
- Trap: the shell blocks `sleep N` chains; poll with `timeout .. bash -c 'until grep -q ...'`. A foreground command over 120 s is moved to the background; use `timeout 100` polls.

## Round 010/013 (T6)
- `rounds/013/toolsmith_verify.py`, `toolsmith_control.py`. n = 9 depth-6 candidates 5e4-6e4 nodes, ~5 ms/node. n = 7 depth-6 control 16/16, depth 5 0/4.

## Round 015 (T6)
- `rounds/015/toolsmith_cords.py N L MINREL PER FIRST LAST [--plan]` (env BOTH, MONO, HIGH, SHORT): n = 8 control with cords; 2 members found at depth 6, not 5 (E-087). Raw visitor members all have sum relations. Cord members come only from LNAs with a 3 early (indices 4, 9-13). Recorded path length is for the labelled algebra: "not at L-1" must be measured.
- Costs: n = 8 depth-6 walk 130 s alone, 200-320 s with 6 in parallel (4+ cores fine); search 250 s.

## Round 017
- Patched `arrowPaths.reduceAgainstPivots` to a full reduction (E-085); test `test_reduce_against_pivots_is_a_normal_form` in tests/test_procedure.py (abstract half fails on old code; real-ideal half passes on both -- the E-085 pair is only checked by `theorist_step7.py`). Patch changes no cord-member count (LNA 4 L = 5: 300, same as before).
- MONO=1 plan n = 8 L = 6 over LNAs 4, 9-13: 0 monomial cord members (`rounds/017/toolsmith_mono_n8_L6_lna*.txt`). Believe (heuristic, untested): a cord needs a two-path (sum) relation in the procedure, so the n = 8 control cannot be monomial and matches the monomial n = 9 candidates only partly.
- Not done: other test files importing procedure/arrowPaths (test_deeper_probing, test_fingerprint, test_invariants, test_lna_moves, test_reflections, test_relation_algebra); LNAs 16-428 never walked for cords.

## Next
- L = 5 MONO `--plan` over LNAs 16-428 in shards; record the producing relation of each cord to test the heuristic.
- Add a unit test on the actual E-085 algebra if the experimentalist's replay supplies a small one.
