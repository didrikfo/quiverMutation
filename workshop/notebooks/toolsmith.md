# Toolsmith notebook

## Round 003
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger keyed by core word; `--jobs`, `--plan`, `--summary`, `--budget-hours`). `45` at 13 pinned in tests/test_orbits_task.py.
- Left out on purpose: a key prefilter that skips walks (E-058: key class is not orbit).
- Ledger is not safe under two concurrent processes on one file; use `--jobs`.

## Round 006 (T3/T8)
- `workshop/rounds/006/toolsmith_orbitclass.py N`: orbit, orbit+mirror, Coxeter-key partitions. Over the 139 placed cores (`--max-word 4`) orbit+mirror refines key everywhere; the 9/10 key-coarser cores are parity classes (E-064, E-070).
- Cost: n = 12 48 s, 13 115 s, 14 598 s, 15 ~11 min, 16 ~32 min (4 slices). Do not run two jobs at once when sizing.
- Trap: the shell blocks `sleep N` chains; poll with `timeout .. bash -c 'until grep -q rc= file; ...'`.

## Round 010 (T6)
- `workshop/rounds/010/toolsmith_verify.py` (`--list`, `--cand I,J`, `--budget-hours H`; exit 2, checked between candidates only). Candidate indices depend on K (K=4: 16 survivors at n = 9, maxdiag -1); name K with the index. Run from the repository root.
- n = 9 depth 5 ~65 s per candidate, depth 6 ~5.5x that.

## Round 013 (T6)
- `workshop/rounds/013/toolsmith_verify.py` adds `nodes` (undeduped visits, both directions) and `distinct` (fingerprint.canonicalKey) to each reached line; `toolsmith_control.py N WALK L PER LO HI` is the depth-L positive control (SHORT=1 for L-1).
- Believe: n = 9 depth-6 candidates are 5.0e4-6.3e4 nodes (4.4e3-7.1e3 distinct), ~5 ms/node; n = 7 depth-6 controls (L = 6 members) find the source 16/16, 3.5e3-1.7e4 nodes; depth 5 finds 0/4. So the E-076 negatives are full searches, but only 4 of 16 measured and the control members are cheap ones (<= 2 relations, mostly hereditary).
- Not known: whether a depth-6 ball can meet a class with no hereditary member at n = 9 (no such control; n <= 8 has none). Depth-7 size: d4->d6 node growth 27x vs time ratio 5.5x per level disagree; time a run.
- Trap: a control job and shards at once on 4 cores inflates seconds, not node counts. `qr.reachedQuipuAlgebras` at n = 7 depth 6 costs 5-50 s per LNA.
- Next: n = 9 non-hereditary-source control (overnight proposal, size one LNA first); still open from 006: n = 12 positive control, `--max-word 5` at n = 14 sizing; round 011's request of A5 reachability (T5, guarded walk from an LNA at n <= 9 to a `tiltingPlus` failure) is the next toolsmith question if nobody has it.
