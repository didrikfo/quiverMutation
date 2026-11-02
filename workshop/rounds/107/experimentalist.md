# T7 / H-010: found no lowering control for a gap-2 pair at k <= 4, margin <= 2 (5 non-vacuous gap-2 patterns, none below base; 4 vacuous), so k = 6 on (8:4)(10:4) was not run -- the shape stays uninformative (no power)

## Response to referee

1. Retitle / reword Claim: done. Title and first Claim sentence now say "found none at k <= 4, margin <= 2", not "does not exist"; the Scope line carries the referee's narrowed wording.
2. Count of non-vacuous gap-2 cells: done, recounted from `experimentalist_chain_out.txt` (no re-run). The 18 patterns listed include 3 that returned None (not evaluated) and 4 run-of-three controls. Of the 11 other evaluated patterns, 4 are start-gap-2 chains that reach no LNA (vacuous: (8,3)(10,3)(12,3), (8,4)(10,4)(12,4), (8,4)(10,4)(12,5), (8,5)(10,5)(12,5)), 5 are start-gap-2 chains that reach an LNA (non-vacuous), and 2 are near-runs with a gap-1 start pair ((8,3)(9,3)(11,3), (8,4)(9,4)(11,4)). I count by pattern, the referee by cell (13 of 18 cells vacuous; the same pattern re-run at another k/margin is a separate cell); the 5 agrees. Correction to "all base -> base": 4 of the 5 are base -> base; (8,4)(10,4)(11,5) goes base 3 -> min 4, i.e. above the base. None is below the base.
3. No conjecture logged: agreed. The Theorist's closing question is withdrawn as a conjecture; there is no mechanism, only an uninformative search.
Also accepted: the 5 non-vacuous patterns are not the stepped chain with the required step count, so the result is "uninformative", not "negative". No re-run was needed; nothing in the numbers changed.


author: experimentalist · round: 107 · kind: negative
thread: T7 · bears on: H-010, F-022, E-230, E-241, E-260, E-271
scope: At k <= 4, margin <= 2, no start-gap-2 chain reached an LNA with overlap below its base; the gap-2 patterns are 4 vacuous (reach no LNA) and 5 non-vacuous (4 base -> base, 1 base -> above base), so no gap-2 control was obtained. Nothing is said about k >= 5. Search: 18 relation patterns (3 returned None, 4 are run-of-three controls) of 3-4 relations (consecutive runs and start-gap-2 chains, lengths 3-5, bases 1-4), k = 3 (margin 1 and 2) and k = 4 (margin 1, four chains only); `lnaMoves.localMutationSequences`, capped by k and margin; "reached 0" cells are vacuous (no LNA inside the window). k = 6 and the (8, 10, 12) chain at k >= 5 were NOT run (cost, below). Capped sample, no proof.
power: this search does say YES where H-010's "run of three" lowering is known: (8,3)(9,3)(10,3) base 2 -> 0 (4 of 6 reached LNAs lowered, k=3 margin 2), (8,4)(9,4)(10,4) base 3 -> 2, (8,5)(9,5)(10,5) base 4 -> 3, (8,3)(9,3)(10,3)(12,3) -> 0. On every gap-2 chain it says nothing: the 4 chains that reach an LNA at all stay at their base. So the harness has power; the gap-2 shapes have no positive control.

## Claim

I found no control of the kind E-271 lacked, at k <= 4 and margin <= 2. For a start-gap-2 pair, (8:4)(10:4) and (8:3)(10:3), no third relation tested (stepped chain with start 12 and length 3, 4 or 5, or 11/12-start variants) gives a window in which an LNA with overlap below the base is reached: at k = 3 and 4 the chains either reach 0 LNAs (13 cells, vacuous), or reach 1-6 LNAs, none below the base overlap (5 non-vacuous patterns: 4 at the base, 1 above it). The only lowerings I find are the consecutive-start runs, which lower the overlap by 1 (from 3 to 2, from 4 to 3) or to 0 for length 3, equal to the E-230 / E-241 picture. Hence a two-bystander cell on (8:4)(10:4) at k = 6 would again be a test that cannot fail, and I did not spend it. This is a statement about k <= 4 and margin <= 2; a lowering may exist at k >= 5 (the stepped chain needs more steps than the E-271 window gave, which is why it is vacuous). It refutes nothing in H-010; a lowering of a gap-2 chain at k = 5-6 would not contradict it either unless it is interior-only.

## Evidence

Full console output: `workshop/rounds/107/experimentalist_chain_out.txt`. Highlights (base -> min over reached LNAs):

| pattern | k, margin | reached | base -> min |
|---|---|---|---|
| (8,3)(9,3)(10,3) | 3, 2 | 6 | 2 -> 0 (control, YES) |
| (8,4)(9,4)(10,4) | 3, 1 | 4 | 3 -> 2 (YES) |
| (8,5)(9,5)(10,5) | 3, 1 | 4 | 4 -> 3 (YES) |
| (8,3)(9,3)(10,3)(12,3) | 3, 1 | 4 | 2 -> 0 (YES) |
| (8,4)(10,4)(12,3) | 3,2 / 3,1 / 4,1 | 2 / 2 / 6 | 2 -> 2 in all three |
| (8,4)(10,4)(11,3 or 4 or 5), (8,3)(10,3)(11,3), (8,3)(9,3)(11,3), (8,4)(9,4)(11,4), (8,4)(10,4)(12,4)(13,4) | 3, 1 | 1-2 | base -> base |
| (8,4)(10,4)(12,4), (8,4)(10,4)(12,5), (8,5)(10,5)(12,5), (8,3)(10,3)(12,3) | 3,1-2 and 4,1 | 0 | vacuous (4-5 steps do not complete it) |

Reading: the stepped chain is the shape whose lowering would need the most steps; at k = 4 it still reaches no LNA, and each added step costs about 2-3 times more (4-step chain cells 3-10 min at margin 1, versus 0.5-1.5 min at k = 3). k = 6 for this shape is hours: an `OVERNIGHT.md` proposal, not a run. Near-gap-1 variants ((8,4)(10,4)(11,x)) look like runs of three but reach only the base.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/107/experimentalist_chain.py 3 2 "8,3:9,3:10,3;8,4:10,4:12,4;8,3:10,3:12,3;8,4:10,4:12,3" --jobs 4   # 78 s
timeout 10m .venv/bin/python workshop/rounds/107/experimentalist_chain.py 3 1 "8,4:9,4:10,4;8,5:9,5:10,5;8,3:9,3:10,3:12,3;8,4:10,4:11,4;8,4:10,4:12,4;8,4:10,4:12,3" --jobs 4   # about 2 min
timeout 10m .venv/bin/python workshop/rounds/107/experimentalist_chain.py 4 1 "8,4:10,4:12,4;8,3:10,3:12,3;8,4:10,4:12,3;8,4:10,4:12,5" --jobs 4   # 9 m 47 s
```

## Prior record

E-271 / `rounds/099/experimentalist.md` Next: "find a lowering control for a start-gap-2 pair ... a stepped chain 8, 10, 12". This round did that search and found none within k <= 4; H-010's run-of-three controls (E-230, E-241, E-260) are reproduced in kind for longer relations (3 -> 2 at length 4, 4 -> 3 at length 5), which I did not find recorded (grep "8,4:9,4" and "9,5" in research/: none), but it is a minor extension, not a new law.

## Code changed

New script only: `workshop/rounds/107/experimentalist_chain.py` (arbitrary relation list, wrapper over `lnaMoves.localMutationSequences`). No library change; no tests run.

## Next

Overnight proposal: `experimentalist_chain.py 5 1 "8,4:10,4:12,3"` and `... 5 1 "8,4:10,4:12,4"` (each likely hours; k = 5 margin 1), to see whether any gap-2 chain lowers; only if one does is a two-bystander cell at k = 6 on (8:4)(10:4) worth running. Theorist: is there a reason a gap-2 chain should be unable to lower (overlap 2 between neighbours leaves no shared-start relation to drop)?
