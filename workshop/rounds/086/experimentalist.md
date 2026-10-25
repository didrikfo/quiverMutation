# Pair (8:3)(9:3), margin 2, k = 4: 5 of 5 sampled m = 2 two-bystander cells (each bystander sharing <= 1 arrow with the pair) do not lower overlap; 3 of 3 completed run-of-three controls do

author: experimentalist · round: 086 · kind: result

## Response to referee

1. Joint run of three (required). Done: `workshop/rounds/086/experimentalist_jointrun.py`, output `experimentalist_jointrun.txt`. A bystander (t,m) covers arrows t..t+m-1; pair covers 8..11. Result per cell, stated plainly: every one of the five cells DOES contain a run of three in bystanders plus pair jointly, because the pair is itself a run of four and every bystander touches or overlaps it (longest joint run 6 or 7). Under the stricter reading (the bystanders' own arrows alone, pair excluded): (11,2)(12,2) YES (arrows 11,12,13), (6,2)(7,2) YES (6,7,8), (7,2)(11,2) no (7,8 | 11,12), (6,2)(11,2) no (6,7 | 11,12), (7,2)(12,2) no (7,8 | 12,13). So 2 of 5 cells (both also the margin-3 cells) have the two bystanders forming a run of three with each other, 3 of 5 do not. The claim is narrowed accordingly: the 5 cells test "no single bystander is a run of three", NOT "no run of three anywhere in the joint cover"; the two overlapping-bystander cells are the E-230 "jointly forming a run" case and they too stay at 2. Controls: all three have a bystander alone covering a run of three or four.
2. Wording (required). Done: "survives" replaced by "is not contradicted by"; claim, title and Scope line use the referee's narrowed wording.
3. 23/30/40 reached LNAs (optional). Not examined; the reached counts stay unexplained.
4. Not done: the controls' least lowering k, k = 5, and the (10,4)(5,2) timeout stay as in Next.

thread: T7 · bears on: H-010, E-230, E-152
scope: For pair (8:3)(9:3), margin 2, k = 4: 5 of 5 sampled m = 2 two-bystander cells (of 18 in E-230's family) with each bystander sharing <= 1 arrow with the pair do not lower overlap; 3 of 3 completed run-of-three controls do. No inference about other cells, k >= 5 or other pair shapes. Detail: pair (8:3)(9:3) only; 4 of the 18 E-230 bystander cells (m = 2 only) at k = 4 margin 2, 2 of them also at k = 4 margin 3 and k = 3 margin 3; 3 controls at k = 4. Capped sample, no proof. Not covered: k >= 5, margin > 3, m = 3 or 4 for the non-control bystanders at k = 4, other pair shapes.

## Claim

With sequences of up to 4 mutations at vertices within 2 of the pattern (`lnaMoves.localMutationSequences`, as E-230), none of five two-bystander cells with every bystander sharing <= 1 arrow with the pair (in 2 of the 5 the bystanders form a run of three with each other; see Response) reaches an LNA of maximum overlap below 2. The same holds at margin 3 (k = 3 and k = 4) for two cells. Controls in which one bystander is a run of three (or a length-4 bystander sharing two arrows) lower the overlap to 0 at k = 4 in 3 of 3 completed cells. So H-010's "only a run of three lowers it" is not contradicted by these cells at k = 4, on the m = 2 cells only. A cell with a bystander of m >= 3 that shares <= 1 arrow was not run at k = 4 (cost). A placement reaching overlap < 2 without a run of three would refute the pattern.

## Evidence

Columns as E-230: base = max overlap of the embedded LNA, reached = LNAs reached, min = least max overlap reached, lowered = reached LNAs with overlap below base.

| steps | margin | bystanders | base | reached | min | lowered | secs |
|---|---|---|---|---|---|---|---|
| 4 | 2 | (11,2)(12,2) | 2 | 23 | 2 | 0 | 175 |
| 4 | 2 | (6,2)(7,2) | 2 | 23 | 2 | 0 | 311 |
| 4 | 2 | (7,2)(11,2) | 2 | 30 | 2 | 0 | 335 |
| 4 | 2 | (6,2)(11,2) | 2 | 40 | 2 | 0 | 505 |
| 4 | 2 | (7,2)(12,2) | 2 | 40 | 2 | 0 | 507 |
| 3 | 3 | (11,2)(12,2) | 2 | 11 | 2 | 0 | 75 |
| 3 | 3 | (6,2)(7,2) | 2 | 11 | 2 | 0 | 73 |
| 4 | 3 | (11,2)(12,2) | 2 | 24 | 2 | 0 | 436 |
| 4 | 3 | (6,2)(7,2) | 2 | 24 | 2 | 0 | 404 |
| 4 | 2 | control (7,3)(12,2) | 2 | 28 | 0 | 15 | 355 |
| 4 | 2 | control (10,3)(5,2) | 2 | 35 | 0 | 18 | 533 |
| 4 | 2 | control (6,4)(12,2) | 2 | 27 | 0 | 6 | 511 |
| 4 | 2 | control (10,4)(5,2) | - | - | - | - | no result in 10 min (timeout, not a verdict) |

Notes. (a) Margin 3 reaches the same 11 LNAs as margin 2 at k = 3, and 24 against 23 at k = 4: widening the window to 3 adds one reached LNA, none lower. (b) Reached counts grow only slowly with k (E-230 k = 3 sibling cells: 11, 15, 23, 11, 28, 23 for the first six m = 2 cells), the cost grows about 4.5x per extra step (175 s against 39 s for (11,2)(12,2)). (c) The run-of-three controls need no more than 4 steps to hit 0 and are not shown to need all 4; I did not record the least k at which they lower (E-230 had them at 3). (d) Test at the cell (6,4)(12,2): bystander (6:4) covers arrows 6..9, sharing arrows 8, 9 with the pair, so it counts as a run of three, matching E-152's `(6:4)`.

## Reproduction

All from the repository root; script `workshop/rounds/086/experimentalist_k4.py` imports `workshop/rounds/081/theorist_twoby.py`. Outputs are `experimentalist_k4_a.txt`, `_k4_ctrl.txt`, `_k4_ctrl4.txt`, `_m3k3.txt`, `_m3k4.txt` in the same directory.

```
timeout 10m .venv/bin/python workshop/rounds/086/experimentalist_k4.py 4 2 "6,2:7,2;7,2:11,2;6,2:11,2;7,2:12,2" --jobs 4   # 8.5 min
timeout 10m .venv/bin/python workshop/rounds/086/experimentalist_k4.py 4 2 "11,2:12,2" --jobs 1                              # 3 min
timeout 10m .venv/bin/python workshop/rounds/086/experimentalist_k4.py 4 2 "7,3:12,2;10,3:5,2" --jobs 2                       # 9 min
timeout 10m .venv/bin/python workshop/rounds/086/experimentalist_k4.py 4 3 "11,2:12,2;6,2:7,2" --jobs 2                       # 7.3 min
timeout 10m .venv/bin/python workshop/rounds/086/experimentalist_k4.py 3 3 "11,2:12,2;6,2:7,2" --jobs 2                       # 1.3 min
timeout 10m .venv/bin/python workshop/rounds/086/experimentalist_k4.py 4 2 "10,4:5,2;6,4:12,2" --jobs 2                       # (6,4) 8.5 min, (10,4) timed out
```

## Prior record

E-230 did k = 3, margin 2, 14 real cells, controls 6 of 6; E-152 did one bystander at k <= 2 fully, k = 3 partially, and lists `(6:4)`, `(10:4)` as lowering by k = 3 (`(10:4)`, `(6:4)` not at k = 2). E-182 states the single-mutation rule that explains why: a slide needs r_v = 2 and an uncovered arrow. Nothing in `research/` goes to k = 4 with two bystanders or to margin 3 (grep "k = 4", "margin 3", "margin > 2": only the "not covered" lists). This is an extension of E-230 along k, not a new kind of statement.

## Code changed

New script only (`workshop/rounds/086/experimentalist_k4.py`, a thin wrapper over the E-230 `run`). No library changes, no tests run.

## Next

Untested and what it costs: k = 5 at about 4.5x more per cell (45 min for the cheapest m = 2 cell; OVERNIGHT, one cell per core); m = 3, 4 non-control bystanders at k = 4 (E-230 k = 3 cells took 150-340 s, so 11-25 min each, over the limit here); margin 4 and above (margin 3 behaved like margin 2 so low priority); other pair shapes at k = 4 (E-152 sampled three at k = 3 only); the 12 other m = 2 E-230 cells at k = 4 (about 4 h of core time). Proposal for OVERNIGHT.md: the 12 remaining m = 2 cells at k = 4 plus the (10,4)(5,2) control, 4 jobs, `timeout` removed. Theorist: a proof of the E-187 type that a lowering sequence needs a run of three would make these cells unnecessary.
