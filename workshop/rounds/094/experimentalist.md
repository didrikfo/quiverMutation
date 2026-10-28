# Pair (8:3)(9:3), k = 5, margin 1 only: in 5 of 5 sampled m = 2 two-bystander cells no LNA with overlap below 2 is reached inside the margin-1 window; 2 of 2 completed run-of-three controls reach 0

author: experimentalist · round: 094 · kind: result
thread: T7 · bears on: H-010, E-230, E-241, E-152
scope: Pair (8:3)(9:3) only; k = 5 mutations; MARGIN 1 (vertices within 1 of the pattern), NOT margin 2 as in E-241 (margin 2 at k = 5 does not fit 10 min). 5 of the 18 E-230 m = 2 cells and 2 of 3 controls. Full enumeration inside the window (no sample cap, no depth cap other than k = 5). Parallel walk checked against the serial walk at k = 5 margin 1 for one cell only ((11,2)(12,2): 31 = 31); the other six cells rely on the k5par reimplementation. In 3 of 5 cells the bystanders alone are not a run of three. Not covered: margin >= 2 at k = 5, other 13 cells, m >= 3 non-control bystanders, other pair shapes, controls (6,4)(12,2) (timeout) and (10,4)(5,2).

## Response to referee

1. Serial check, DONE. `experimentalist_k5.py 5 1 "11,2:12,2" --jobs 1` (the library `localMutationSequences` walk): reached 31, min 2, lowered 0, 213 s. The parallel run gave 31, min 2. They agree at k = 5 margin 1 for this cell. The other six cells are still checked only against k5par (k = 3, 4 margin 2 plus this cell).
2. Run of three, DONE. The Claim now says that only 2 of 5 cells have bystanders alone forming a run of three; for the other 3 the test is only that no single bystander is one. The title and Scope are narrowed to "inside the margin-1 window" and "k5par reimplementation".
3. Margin 1 at k = 4, DONE. `experimentalist_k5par.py 4 1 "11,2:12,2" --jobs 3`: reached 21, min 2, lowered 0, 35 s. Against E-241 margin 2: 23 reached. So margin 1 loses 2 LNAs of 23 at k = 4 in this cell (none lower). It is not shown that margin 1 loses nothing at k = 5; the lost LNAs at k = 4 have overlap 2.
4. Wording: the Claim opens with the window and the reimplemented walk; "H-010 is not contradicted" kept.

New table row (k = 4, margin 1): see Evidence. Note (c) is superseded by item 3.

## Claim

Within the margin-1 window, and with the k5par reimplementation of the walk, with sequences of up to 5 mutations at vertices within 1 of the pattern's arrow span, none of five m = 2 two-bystander cells (each bystander sharing <= 1 arrow with the pair; in only 2 of the 5, (11,2)(12,2) and (6,2)(7,2), do the bystanders alone form a run of three, so for the other 3 the test is only that no single bystander is a run of three) reaches an LNA of maximum overlap below 2. Controls with a bystander that is a run of three, (7:3)(12:2) and (10:3)(5:2), lower it to 0 (21 of 37 and 25 of 46 reached LNAs lower). So H-010's "only a run of three lowers the overlap" is not contradicted at k = 5, margin 1, on these cells. This is a narrower window than E-241 (margin 2): the margin-1 k = 5 result says nothing about sequences that need a vertex two away. A placement reaching overlap < 2 with no run of three would refute the pattern; none found. The k = 5 margin-2 cells did not complete (see Evidence), which is a cost verdict, not a result.

## Evidence

Columns as E-241. `secs` = wall clock of `experimentalist_k5par.py` with the stated jobs (two cells run concurrently at 2 jobs each, so times are not comparable across rows).

| steps | margin | bystanders | base | reached | min | lowered | secs (jobs) |
|---|---|---|---|---|---|---|---|
| 5 | 1 | (11,2)(12,2) | 2 | 31 | 2 | 0 | 156 (4) |
| 5 | 1 | (6,2)(7,2) | 2 | 31 | 2 | 0 | 256 (2) |
| 5 | 1 | (7,2)(11,2) | 2 | 38 | 2 | 0 | 223 (2) |
| 5 | 1 | (6,2)(11,2) | 2 | 57 | 2 | 0 | 567 (2) |
| 5 | 1 | (7,2)(12,2) | 2 | 57 | 2 | 0 | 597 (2) |
| 5 | 1 | control (7,3)(12,2) | 2 | 37 | 0 | 21 | 590 (2) |
| 5 | 1 | control (10,3)(5,2) | 2 | 46 | 0 | 25 | 527 (4) |
| 5 | 1 | control (6,4)(12,2) | - | - | - | - | timeout at 13 of 14 first moves, 43 reached so far, no verdict |
| 5 | 2 | (11,2)(12,2) | - | - | - | - | serial script timed out at 600 s (est. 790 s from E-241 x 4.5) |
| 4 | 1 | (11,2)(12,2) | 2 | 21 | 2 | 0 | 35 (3), margin 1 vs E-241's 23 at margin 2 |
| 5 | 1 | (11,2)(12,2) serial | 2 | 31 | 2 | 0 | 213 (1), library walk, equals k5par |
| 4 | 2 | (11,2)(12,2) | 2 | 23 | 2 | 0 | 207 (4), parallel script reproduces E-241's 23 (serial 175 s) |

Plan (`--plan`, from E-241 k = 4 seconds x 4.5): margin-2 k = 5 estimates are 790 s (11,2)(12,2) to 2400 s (10,3)(5,2), i.e. none fits one 10-minute command; the parallel split gave no speed-up (207 s against 175 s at k = 4) because each worker loses the shared `seen` table. Margin 1 shrinks the allowed vertex set by 2 and brings every cell except (6,4)(12,2) under 10 min.
Notes. (a) The k = 5, margin 1 reached counts (31..57) exceed E-241's k = 4, margin 2 counts (23..40): more steps matter more than the margin here. (b) 2 of the 5 cells ((11,2)(12,2), (6,2)(7,2)) have the bystanders forming a run of three with each other (E-241 Response); they too stay at 2. (c) I did not run margin 1 at k = 4, so I cannot say whether margin 1 loses LNAs that margin 2 finds at k = 5; at k = 3 and 4 E-241 reports margin 3 adds one LNA and none lower. (d) Controls not shown to need 5 steps; least lowering k not recorded (E-230: 3).

## Reproduction

From the repository root. `experimentalist_k5par.py` splits the library's `lnaMoves.localMutationSequences` search by first move (same pruning inside each branch; union equals the serial result, checked at k = 3 margin 2: 11 reached, and k = 4 margin 2: 23 reached).

```
.venv/bin/python workshop/rounds/094/experimentalist_k5.py 5 2 "11,2:12,2;6,2:7,2;7,2:11,2;6,2:11,2;7,2:12,2;7,3:12,2;10,3:5,2;6,4:12,2" --plan   # instant
timeout 10m .venv/bin/python workshop/rounds/094/experimentalist_k5par.py 5 1 "11,2:12,2" --jobs 4    # 156 s
timeout 10m .venv/bin/python workshop/rounds/094/experimentalist_k5par.py 5 1 "6,2:7,2" --jobs 2      # 256 s (run with the next one)
timeout 10m .venv/bin/python workshop/rounds/094/experimentalist_k5par.py 5 1 "7,3:12,2" --jobs 2     # 590 s
timeout 10m .venv/bin/python workshop/rounds/094/experimentalist_k5par.py 5 1 "7,2:11,2" --jobs 2     # 223 s
timeout 10m .venv/bin/python workshop/rounds/094/experimentalist_k5par.py 5 1 "6,2:11,2" --jobs 2     # 567 s
timeout 10m .venv/bin/python workshop/rounds/094/experimentalist_k5par.py 5 1 "7,2:12,2" --jobs 2     # 597 s
timeout 10m .venv/bin/python workshop/rounds/094/experimentalist_k5par.py 5 1 "10,3:5,2" --jobs 4     # 527 s
timeout 10m .venv/bin/python workshop/rounds/094/experimentalist_k5.py 5 2 "11,2:12,2" --jobs 1       # times out (exit 124)
```
Outputs: `workshop/rounds/094/experimentalist_k5_m1_*.txt`, `experimentalist_k5_a.txt` (the timeout), all under 1 KB.

## Prior record

E-241 (k = 4, margin 2, 5 cells, 3 controls), E-230 (k = 3), E-152/E-153 (one bystander), E-182/E-187 (single-mutation rule). grep of `research/EXPERIMENTS.md` finds no k = 5 result for H-010. This extends E-241 along k at a narrower margin; it is not a new kind of statement and not a proof.

## Code changed

New scripts only: `experimentalist_k5.py` (copy of 086 k4 plus `--plan`), `experimentalist_k5par.py` (first-move-split reimplementation of the `localMutationSequences` walk using the library's private helpers `_quiet`, `_copy`, `_stateKey`, `asRelLengths`). No library changes; no tests run.

## Next

OVERNIGHT proposal: k = 5 at margin 2 for the 5 non-control cells plus 3 controls, serial `localMutationSequences` (shared `seen` table), 790-2400 s each by the x4.5 estimate, one cell per core, 4 cores, about 1 h wall; and (6,4)(12,2), (10,4)(5,2) at k = 5 margin 1. Toolsmith: a faster `localMutationSequences` (state-key caching or a shared seen table across processes) would make k = 6 reachable. Theorist: a proof that a lowering sequence needs a run of three would retire this family.
