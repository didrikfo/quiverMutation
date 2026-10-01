# Toolsmith notebook

## Round 003
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger keyed by core word; `--jobs`, `--plan`, `--summary`, `--budget-hours`). `45` at 13 pinned in tests/test_orbits_task.py.
- Left out on purpose: a key prefilter that skips walks (E-058: key class is not orbit).
- Ledger is not safe under two concurrent processes on one file; use `--jobs`.

## Round 006 (T3/T8)
- `workshop/rounds/006/toolsmith_orbitclass.py N` compares orbit, orbit+mirror and Coxeter-key partitions. Over the 139 placed cores
  (`--max-word 4`) orbit+mirror refines key everywhere; the 9/10 key-coarser cores are parity classes (E-064, E-070).
- Cost: n = 12 48 s, 13 115 s, 14 598 s, 15 ~11 min, 16 ~32 min (4 slices). A shared machine slows slices; do not run two jobs at once when sizing.
- Trap: the shell blocks `sleep N` chains; poll with `timeout .. bash -c 'until grep -q rc= file; ...'`.

## Round 010 (T6)
- The H-017 candidate script is `workshop/rounds/004/maverick_verify.py` (round 009 had no copy). Now
  `workshop/rounds/010/toolsmith_verify.py` with `--list`, `--cand I,J`, `--budget-hours H` (exit 2, checked between candidates only).
- Believe: candidate indices depend on K (K=1: 4, K=4: 16, K=100: 160 survivors at n = 9, maxdiag -1); name K with the index.
- n = 9 depth 5, K = 1: reached [] x4 (261 s); depth 6 cand 2: reached [] in 434 s (ratio 5.5 to depth 5). One candidate per 10-minute shard fits, margin 1.4x.
- Run from the repository root: the script needs `families.py` on the path (it adds cwd).
- Not done: depth 6 for the other 15 K = 4 candidates (12 untimed at depth 5); an in-candidate budget hook; any test (script is under workshop/, untested).
- Next: if the chair has slots, run the 16 shards and log which exceed the cap; then propose depth 7 to OVERNIGHT.md. Also still open from 006: n = 12 positive control for H-017, `--max-word 5` at n = 14 sizing.
