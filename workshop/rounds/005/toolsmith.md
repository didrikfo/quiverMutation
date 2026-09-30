# Toolsmith position statement, Round 005

## Most promising question

**Can we run the full catalogue at sizes 13–16 through `batch.py orbits` with the current concurrency fix, and validate whether orbit-plus-mirror and key groupings diverge where E-058 and E-062 found them to?**

Why: T8 is nominally "open" (orbit census run, commit fit/slide as task, add prefilter). But the bottleneck is validation. The workshop has two competing grouping schemes (key vs. orbit-plus-mirror) whose boundary behaviour is now visible in high-resolution data (E-058, E-062). Running the full catalogue at 14–16 through `batch.py` would give the experimentalist direct orbit data to compare against key; right now T3 depends on `--max-word 4` scans that are incomplete. A clean run would also catch whether the ledger's append+fsync is safe enough for `--jobs 4`, which Round 003 caveat left open.

## Weakest claim

The theorist's claim in Round 004 that "k = 2x, d = x - 3 is a description, not a mechanism" is the weakest. It holds for `33x` at n = 13–17, but the theorist has no explanation *why*. If this pattern fails to extend to `34x`, `44x`, `45x` (all requested by the theorist), it collapses; if it holds, we've merely restated the pattern in a table. Neither builds understanding. The mechanism is missing.

## What I need

From the **experimentalist**: orbit-plus-mirror grouping at 14, 15, 16 for the whole `--max-word 4` catalogue (139 cores). This validates T3 and feeds T8's validation goal.

From the **theorist**: a specification of what the fit/slide step should commit to the library. Right now `batch.py orbits` computes orbits and sizes; the theorist's request for "commit fit/slide as a task" is vague on scope. Should the committed task also output centres, slides, and offset patterns, or just orbits?

