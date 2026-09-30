# Toolsmith notebook

## Round 003

- Task T8: orbit report of E-052 as a `batch.py` task. Done; see
  `workshop/rounds/003/toolsmith.md`.
- Choice: put `orbitCensus` in `batch.py` (needs `_rowFor`, which lives there)
  rather than the library; the census script's JSONL/shard resume is replaced by
  `jobs.Ledger` keyed by core word, so `--jobs`, `--plan`, `--summary`,
  `--budget-hours` all work for free.
- The old script's logic is sound as far as I can tell; the only oddity is its
  docstring path (says rounds/001) and that its output is not the ledger format.
- Timings: `45` at 13 is 8 s, so the pin is a plain non-slow test; no smaller
  companion needed beyond the n = 9 resume test.
- Left out on purpose: a Coxeter-key prefilter. E-058 says key class is not
  orbit; the only safe prefilter is orbit-plus-mirror, and nothing needs it yet.
- Caveat noticed: the ledger is not safe under two concurrent processes on the
  same file (append with fsync, but each process computes its own "done" set at
  start); use `--jobs`, not several shells, on one ledger.
- Did not run the whole catalogue at 13 through the new task; only three cores.
  Round 002 has the 139-core result; a full rerun through the task would confirm
  the two agree byte for byte in the fields.
