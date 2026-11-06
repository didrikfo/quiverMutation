# No gap-2 pair lowering control: three non-vacuous pair cells (k = 4: two, k = 5: one), none lowered; k = 5, 6 otherwise capped at 170 s, (8,3)(10,3) hits an engine ValueError at k = 5 (no power)

author: experimentalist · round: 119 · kind: negative
thread: T7 · bears on: H-010, E-288, E-230, E-241, E-260, E-271
scope: start-gap-2 pairs only, margin 1; k = 4: 4 of 5 cells finished, of which 2 non-vacuous; k = 5: 1 non-vacuous cell ((8,4)(10,3), uncapped 223 s); k = 5, 6 otherwise capped at 170 s per cell or ValueError, not measured beyond that; no bystanders; k >= 6 not decided
power: none -- no cell lowers overlap; the harness can say YES (E-288: runs of three, 2 -> 0, 3 -> 2) but no gap-2 control exists, so "stays at 2" cannot fail (no power)
backlog: item 10 (T7 / H-010): converted (overnight proposal below); not settled, not retired

## Response to referee

1. k=5 row for (8,4)(10,3): fixed to CAP 170 s (what `_k5.txt` records); the k=6 column now says what its 3-line file has, and (8,5)(10,4) is "not run". Done. Also re-ran (8,4)(10,3) at k=5 uncapped (`experimentalist_trace.py`, `experimentalist_trace_k5_84_103.txt`): 11 LNAs, min 2, 0 lowered, 223 s. That is a third non-vacuous cell; the title/scope/claim now say so.
2. "4 of 5 finished" now reads "4 of 5 finished, 2 non-vacuous" at k=4 (scope line and Evidence note define vacuous). Vacuous cells: (8,4)(10,4) and (8,3)(10,3).
3. E-288 comparison: E-288 ran three-relation chains (8,10,12 and variants), not two-relation pairs. Its 5 non-vacuous gap-2 patterns all carry a third relation, e.g. (8,4)(10,4)(12,3); the pairs (8,4)(10,3) and (8,5)(10,4) are not among its listed patterns, so they are new as patterns but not as a verdict (still base -> base). New evidence for H-010: zero lowerings, no change of verdict; new content is the bystander-free shape. I did not match E-288's patterns further than its output file `rounds/107/experimentalist_chain_out.txt`.
4. "hours each" removed. Measured: 34-170 s per pair at k=4, 223 s for one k=5 pair; the rest is "unmeasured, > 170 s". The proposal now states that E-288 already declined it, and recommends not approving it again until the theorist question is answered (I give no reason to approve; k=5 for one cell is cheap, k=6 unmeasured).
5. ValueError diagnosed with a traceback (59.5 s): engine limitation by documented design (`arrowPaths.liftPath` cannot read a vertex-sequence relation across parallel arrows, called from `procedure.relationsFrom` during the mutation walk at depth 4); not wrapper misuse. Whether the walk should skip such branches is a library decision, not made here.
6. Theorist question (does H-010's end-only argument exclude lowering for gap-2 starts?): recorded below as the preferred retirement route; for next round, not done by me.
7. Not done: re-running k=6 or the other k=5 pairs (each > 170 s, unmeasured; would not change a verdict that has no power).

## Claim

I found no start-gap-2 pair that lowers overlap, in three non-vacuous cells (k = 4: (8,4)(10,3), (8,5)(10,4); k = 5: (8,4)(10,3)); k = 6 was reached for no cell within a 170 s cap. Pairs (8,l)(10,l') with margin 1 stay at their base overlap whenever they reach an LNA: the two-bystander shape (E-230/E-241/E-260/E-271) still has no control. This does not show that no lowering gap-2 shape exists at k >= 6; k = 5, 6 were undecided, not negative. Refuted by any gap-2 cell that reaches an LNA with overlap below base.

## Evidence

Harness: E-288's `localMutationSequences` via `rounds/107/experimentalist_chain.run`, wrapped in `rounds/119/experimentalist_gap2.py` (per-cell alarm, exceptions recorded). Pairs only, margin 1.

| pair | base | k = 4 | k = 5 | k = 6 |
|---|---|---|---|---|
| (8,4)(10,4) | 2 | 0 LNAs reached (109 s), vacuous | CAP 170 s | CAP 170 s |
| (8,5)(10,5) | 3 | CAP 170 s | CAP 170 s | CAP 170 s |
| (8,3)(10,3) | 1 | 0 LNAs (34 s), vacuous | ValueError (parallel arrows 10->13) | ValueError (8->11) |
| (8,4)(10,3) | 2 | 6 LNAs, min 2, 0 lowered (49 s) | CAP 170 s (uncapped re-run: 11 LNAs, min 2, 0 lowered, 223 s) | not run in the wrapper (k=6 file has 3 lines) |
| (8,5)(10,4) | 3 | 4 LNAs, min 3, 0 lowered (168 s) | not run in the wrapper | not run in the wrapper |

("Not run": the cell was still queued when the 10-minute wrapper ended; the k=5 file has 4 lines, the k=6 file 3. Vacuous = no LNA reached, so the cell cannot show a lowering and has no power even for "stays at base".)
Sizing (measured): one pair at k = 4, margin 1 costs 34-170 s; one k = 5 pair took 223 s uncapped ((8,4)(10,3)); (8,4)(10,4), (8,5)(10,5) at k = 5, 6 exceeded 170 s and are unmeasured beyond that; (8,5)(10,4) was not run at k = 5, 6. ValueError diagnosis (traceback in `experimentalist_trace_k5_83_103.txt`, 59.5 s): the error comes from `lnaMoves.localMutationSequences.walk` -> `mutation.quiverMutationAtVertices` -> `procedure.relationsFrom` -> `arrowPaths.lift`, i.e. the engine re-reads stored vertex-sequence relations as arrow paths after 4 mutation steps, and `liftPath`'s docstring states that this restriction (no parallel arrows along a stored vertex path) is by design; so it is an engine limitation reached by the walk, not wrapper misuse (the wrapper only calls the library entry point), though whether the library should instead skip such a branch is unjudged. A cap is recorded as undecided, never as outside.

## Reproduction

```
for k in 4 5 6; do timeout 10m .venv/bin/python workshop/rounds/119/experimentalist_gap2.py $k 1 170 "8,4:10,4;8,5:10,5;8,3:10,3;8,4:10,3;8,5:10,4" > workshop/rounds/119/experimentalist_gap2_k$k.txt; done   # about 10 min each run in parallel
```

Output files `experimentalist_gap2_k4.txt`, `_k5.txt`, `_k6.txt` (a few lines each).

## Overnight proposal (for OVERNIGHT.md)

E-288 already declined an overnight proposal of this kind ("not approved"). This one adds nothing new except that k = 5 for (8,4)(10,3) is measured at 223 s, so a k = 5 run is minutes, not hours; I propose it not be approved again for k = 6 until the theorist route below is answered. If approved anyway: run `experimentalist_gap2.py 5 1 3000` and `6 1 6000` on (8,4)(10,4), (8,5)(10,4) as separate background jobs (unmeasured, each > 170 s), (the parallel-arrow ValueError only hits (8,3)(10,3); the wrapper records it). A YES on any gap-2 cell is the control; if none, T7 retires with scope "gap-2 pairs, k <= 6, margin 1". A cheaper route: ask the theorist whether a gap-2 pair can be shown unable to lower (e.g. lowering needs a gap-1 start adjacent to a mutated vertex), which would retire the item with a reason.

## Prior record

E-288 (no gap-2 control at k <= 4, no power), E-271 (shape (8:4)(10:4), no control), E-230/E-241/E-260. This repeats E-288's null on pairs without bystanders and adds sizing, one k = 5 cell and the k >= 5 engine error; nothing contradicts H-010.

## Code changed

Added `workshop/rounds/119/experimentalist_gap2.py` (wrapper) and `experimentalist_trace.py` (traceback and uncapped timing); no library change, no tests run.
