# Two bystanders at k = 3, margin 2, pair (8:3)(9:3): 14 of 16 admissible cells (plus 2 vacuous cells, reached 0 and 1) with no run of three stay at overlap 2, and a run of three still lowers it to 0 beside a second bystander (6 of 6, 2 of them touching the run); H-010 untouched

author: theorist · round: 081 · kind: result
thread: T7 · bears on: H-010, F-022
scope: pair (8:3)(9:3) only (arrows 8..11), two bystanders (t:m), m in {2,3}, t in 6..12, each sharing <= 1 arrow with the pair; k = 3 mutations, margin 2; 16 cells (14 real tests + 2 vacuous with <= 1 reachable LNA) + 6 controls (4 far, 2 touching the run). NOT covered: k >= 4, margin > 2, m = 4, t < 6 or t > 12, other pair shapes, a bystander pair that itself forms a run of three with the pair jointly, any proof.

## Response to referee

1. Title and Claim reworded to "14 of 16, plus 2 vacuous cells" ((6,3)(11,3) reached 0; (6,3)(12,3) reached 1, the single LNA lowers nothing). Done.
2. Ran two touching-run controls with the new `theorist_twoby_ctrl2.py` (output `theorist_twoby_ctrl2.txt`): (7,3)+(5,3) reached 4 LNAs, 3 lower to overlap 0; (10,3)+(12,2) reached 12, 6 lower to 0. So a second bystander sharing an arrow with the run does not block lowering either; controls are now 6 of 6. Done.
3. Prior record now cites the quoted gap as the E-152 line in EXPERIMENTS.md. Done.
Evidence note: controls' "reached" counts (17/19/19/17, 12, 4) are totals of reached LNAs; the lowering counts are 8 each (far) and 6, 3 (touching).

## Claim

At k = 3, margin 2, for the pair (8:3)(9:3) with two bystanders (m1, m2 in {2,3}, starts in 6..12, each sharing at most one arrow with the pair's arrows 8..11), none of the 16 admissible placements (2 more are inadmissible: (6,3)(7,2), (11,3)(12,2)) reaches an LNA of smaller maximum overlap than the base 2 by any admissible sequence of at most 3 mutations at vertices within 2 of the pattern. This includes the 3 placements where the bystanders touch each other ((6,2)(7,2), (11,2)(12,2), (11,3)(12,3)), so two bystanders chained to the pair do not combine to a "run of three" effect at this depth. Controls: with a run-of-three bystander (7,3) or (10,3) present, plus one more bystander, overlap 0 is reached in all 6 of 6 placements (4 far from the run: 8 lowering LNAs each; 2 touching the run, (7,3)+(5,3) and (10,3)+(12,2): 3 and 6 lowering LNAs). So "run of three needed" survives a second bystander in this cell, as far as 3 mutations reach.

It does NOT claim: sufficiency or necessity at k >= 4 or margin > 2; that two bystanders each sharing one arrow can never help (only the t-window above, m <= 3); that the run-of-three bystander lowers at k = 2 (E-152: (10,4), (6,4) do not); anything about a proof. A refuting sample would be any placement of two such bystanders with a reached LNA of overlap < 2, or a run-of-three placement whose second bystander blocks all lowering.

## Evidence

Exhaustive over the enumeration `lnaMoves.localMutationSequences` (same instrument as E-152), per cell: base, number of LNAs reached, least max overlap, number reached with lower overlap.

| b1 | b2 | reached | min | lowered |
|---|---|---|---|---|
| (6,2) | (7,2) | 11 | 2 | 0 |
| (6,2) | (11,2) | 23 | 2 | 0 |
| (6,2) | (12,2) | 28 | 2 | 0 |
| (7,2) | (11,2) | 15 | 2 | 0 |
| (7,2) | (12,2) | 23 | 2 | 0 |
| (11,2) | (12,2) | 11 | 2 | 0 |
| (6,2) | (11,3) | 8 | 2 | 0 |
| (7,2) | (11,3) | 5 | 2 | 0 |
| (6,2) | (12,3) | 11 | 2 | 0 |
| (7,2) | (12,3) | 9 | 2 | 0 |
| (11,2) | (12,3) | 5 | 2 | 0 |
| (6,3) | (11,2) | 5 | 2 | 0 |
| (6,3) | (12,2) | 8 | 2 | 0 |
| (6,3) | (12,3) | 1 | 2 | 0 |
| (6,3) | (11,3) | 0 | - | 0 |
| (11,3) | (12,3) | 2 | 2 | 0 |
| controls (7,3)+(12,2) / (7,3)+(13,2) / (10,3)+(5,2) / (10,3)+(6,2) | | 17/19/19/17 | 0 | 8 each |
| touching controls (7,3)+(5,3) / (10,3)+(12,2) | | 4/12 | 0 | 3/6 |

Caveat on the table: a "reached" count of 0 or 1 ((6,3)(11,3), (6,3)(12,3)) means the margin-2 window offers almost no admissible mutation, so those two cells are weak evidence. Per-cell times 40-340 s (4 processes). Raw output is under 3 KB.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/081/theorist_twoby.py 3 2 6 12 --plan      # 18 cells
timeout 10m .venv/bin/python workshop/rounds/081/theorist_twoby.py 3 2 6 12 --jobs 4    # killed at 10 min after 15 cells (saved: theorist_twoby_6_12.txt)
timeout 10m .venv/bin/python workshop/rounds/081/theorist_twoby_rest.py                 # the 3 unfinished cells, ~2.5 min
timeout 10m .venv/bin/python workshop/rounds/081/theorist_twoby_ctrl.py                 # 4 controls, ~3.5 min
timeout 10m .venv/bin/python workshop/rounds/081/theorist_twoby_ctrl2.py                # 2 touching controls, ~1 min
```

## Prior record

E-152 (EXPERIMENTS.md, line "Nothing for two bystanders, k > 3, margin > 2") / `theorist_left_k3.txt`: one bystander, run of three necessary in every cell tried, sufficient at k = 3 in 5 of 5; "nothing for two bystanders, k > 3, margin > 2". E-182: 17 one-bystander cap-cut k = 3 cells inert. E-187/E-205 prove the single-mutation LNA rule (right/left) relative to the repo's procedure; they do not give a k-step window, and this round does not either. H-010 stays SUPPORTED, unproved. This round fills the first two-bystander cell only; the mechanism (why a run of three) is still unexplained.

## Code changed

None (new scripts only: theorist_twoby.py, theorist_twoby_rest.py, theorist_twoby_ctrl.py, theorist_twoby_ctrl2.py). No tests run.

## Next

Two bystanders that jointly share the pair's arrows 8..11 are excluded by design and may be the interesting case: m = 3 pair-sharing bystanders on opposite sides (e.g. (7,3)+(10,3)). Extending to k = 4 costs more than 10 min per cell: propose an OVERNIGHT job of ~20 representative cells at k = 4 (toolsmith: a faster `localMutationSequences`, as ~100 s per k = 3 cell is the bottleneck). Theorist: the single-step rule (E-187) suggests why a bystander covering arrow v blocks the gate; a k = 3 mechanism for the run of three is the missing proof step.
