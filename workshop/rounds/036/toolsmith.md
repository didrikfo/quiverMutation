# Toolsmith round 036

## Most promising question

Why does the n=7 both-die BFS frontier grow exponentially (ratio 2.4–2.6) instead of falling? This is the gate for all downstream: forward search exhausts budget, and reverse search feasibility depends on frontier size.

**Why it matters:** E-130 sized the forward closure: class 1 ratio 2.2→2.5 (29k at level 8, 240 s), class 3 ratio 2.6→1.9 (40k at level 7), 0/44 targets reached in 10 min each. Frontier growth is the bottleneck. If the ratio is structural (mutation algebra property), we need a different strategy (hand-built examples, guessing from invariants, or acceptance of non-reachability). If it's an artifact (bad BFS order, key drops, repeats in the seen set), it's a toolsmith fix.

## Weakest claim the workshop relies on

**"0/44 reached in 240 s per class" assumes all 44 targets have non-None Coxeter keys.** My notebook (round 035) flags: does `_coxeterKeyOrNone` return None on any target? If yes, the BFS skips them silently, and the count is misleading. This is E-124's base case and must be verified before any reverse search.

## What I need from another persona

- **Experimentalist:** diagnostic run on the 44 target states (n=7, both-die classes 1 and 3). Do all have non-None keys? Any orbit or key patterns across the 44?
- **Theorist:** why the frontier ratio stays high: is it a structural property of mutations from both-die roots, or an artifact (e.g., key collisions merging distinct mutation paths)?

---

*Next: Reverse search from the 44 targets, sized with frontier data and key diagnostics. If frontier stays >2, hand-build examples instead.*
