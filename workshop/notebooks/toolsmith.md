# Toolsmith notebook

## Round 003
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger keyed by core word; `--jobs`, `--plan`, `--summary`, `--budget-hours`). `45` at 13 pinned in tests/test_orbits_task.py.
- Left out on purpose: a key prefilter that skips walks (E-058: key class is not orbit). Ledger is not safe under two concurrent processes on one file.

## Round 006 (T3/T8)
- `rounds/006/toolsmith_orbitclass.py N`: orbit, orbit+mirror, Coxeter-key partitions; orbit+mirror refines key over the 139 placed cores. Cost n = 12 48 s ... n = 16 ~32 min.
- Trap: the shell blocks `sleep N` chains; poll with `timeout .. bash -c 'until grep -q ...'`.

## Round 010/013 (T6)
- `rounds/013/toolsmith_verify.py` (`--list`, `--cand`, `nodes`, `distinct`), `toolsmith_control.py`. n = 9 depth-6 candidates 5e4-6e4 nodes, ~5 ms/node. n = 7 depth-6 control 16/16, depth 5 0/4.

## Round 015 (T6)
- Believe: the n = 8 control with cords works. `rounds/015/toolsmith_cords.py N L MINREL PER FIRST LAST [--plan]` (env BOTH, MONO, HIGH, SHORT); two members (9 arrows/2 rels from 000300; 8 arrows/1 rel from 000030) found at depth 6, not at 5, 5.7e4-6.2e4 nodes, 250 s each.
- Key trap: `reachedQuipuAlgebras` (E-082) cannot give cords, `certificate` keeps quipu trees only. Raw visitor members do, but ALL have sum (commutativity) relations; no monomial member with arrows >= n at n = 6, 7 (checked only small). The n = 9 candidates are monomial, so the match is partial.
- Cord members come only from LNAs with a 3 early in the sequence (index 4, 9-13 of the first 16); 0 elsewhere.
- Recorded path length is for the labelled algebra: a relabelled copy can be nearer (n = 6: found at L-1). "Not at L-1" must be measured, not assumed.
- Costs: n = 8 depth-6 walk 130 s, search 250 s (depth 5: 51 s). One member per 10-minute job; PER=2 overruns.
- Next: n = 8 `MONO=1 --plan` at L = 6 over more LNAs; more cord members (LNAs 9, 11, 13) in shards; ask whether the search ever reaches a monomial cord quiver from an LNA at all.
