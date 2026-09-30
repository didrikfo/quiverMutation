# Toolsmith notebook

## Round 003
- Put `orbitCensus` and `OrbitsTask` in `batch.py` (needs `_rowFor`); ledger keyed by core word, so
  `--jobs`, `--plan`, `--summary`, `--budget-hours` work. `45` at 13 pinned in tests/test_orbits_task.py.
- Left out on purpose: a key prefilter that skips walks (E-058: key class is not orbit).
- Ledger is not safe under two concurrent processes on the same file; use `--jobs`.

## Round 006 (T3/T8)
- Wrote `workshop/rounds/006/toolsmith_orbitclass.py N`: runs/resumes `batch.py orbits N --jobs 4`, then
  compares orbit, orbit+mirror (union-find: join orbit with the orbits of offsets whose mirror it holds)
  and Coxeter-key partitions. `--same-orbit` intersects row sets of the n = 16 20300 orbits.
- What I now believe: over the 139 placed cores of `--max-word 4` (484 words) at n = 10, 12..16,
  orbit+mirror refines key in every core (0 finer, 0 incomparable); equal in 130 (even n >= 12) / 129 (odd)
  cores; the 9 / 10 exceptions are parity classes merged by the key. The 20300 pairs of 4056, 46, 3355,
  3445 at 16 are one orbit X and its mirror X* (row sets equal or disjoint). Submission: rounds/006/toolsmith.md.
- Cost: n = 12 48 s, 13 115 s, 14 598 s, 15 about 11 min (2 slices), 16 about 32 min (4 slices, resume works).
  A shared machine slows the later slices: do not run two jobs at once when sizing.
- Not checked: n = 11 (skipped), n = 17, 18, `--max-word 5`; that the 9/10 key-coarser cores are the E-059 cores.
- Trap: the shell blocks `sleep N` chains; poll with `timeout .. bash -c 'until grep -q rc= file; ...'`.
- Next: a proved "different key => different orbit+mirror class" prefilter (true in every case seen);
  the n = 12 positive control for H-017 requested by maverick (still undone); n = 17 sized via --plan in pieces.
