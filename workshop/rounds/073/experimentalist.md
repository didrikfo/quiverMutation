# S-1 at n = 10 -> 9: the E-190 pair is outside the D >= 3 rule's domain; at n = 10 -> 9 all deletions of 90000000 and 7 of 10 of 00000030 land in P^(6)_(1,1), the rest in the hereditary class. Exact classifyLength(10) is sized, not run

author: experimentalist · round: 073 · kind: result
thread: S-1 · bears on: H-018, H-020 (S-1 only)
scope: n = 10 -> 9, one pair (two LNAs, all 10 deletions each); n = 9 labels exact (round-065 classifyLength(9) cache, 1430 rows); plan for n = 10 is a size estimate from 6 timed searches, not a run; image classes compared by label only (not as one n = 9 object)

## Response to referee

1. Retitled as required; "under unrestricted deletion both images lie in one class" dropped from title and claim. The claim now says plainly that "same class for i = 1..7" is a fact about 00000030's own deletions (7 of 10), not evidence for a transport rule. Done.
2. Convention stated (Claim, (b)): `removeVertex` contracts through the deleted vertex. Example: deleting i = 7 from the relation 7..10 of 00000030 gives a length-3 relation on 6..9, i.e. 00000300; E-206 uses the same convention. Done.
3. Loop fixed. Exit 124 (a class step overran `timeout`) now also resumes; the CSV and progress sidecar are resumable after a kill. New loop (tested with stub commands: exit 124 and 2 continue, 0 and 1 stop, six passes cap; and twice with the real `--run` at n = 8, exit 2 both times, unplaced 10 then 8): see Proposal. Worst case is 6 x 9.5 h = 57 h; stated. Done.
4. "~1 h to > 9 h" replaced by "unknown". The 3 h figure is quoted only as the arithmetic of 6 evenly spaced samples, not a bound. Done.
5. Whether 80000000 and 00000300 are joined by orbit+mirror moves at n = 9: next round (not done this sitting). Note both carry the label P^(6)_(1,1) in the exact cache, so they are in one class; what is open is the orbit+mirror question, which the Next section keeps.

Also added (hole 3): 90000000 is degenerate, every deletion of a rad^9 algebra is rad^8 at n = 9, so its column of the table is trivial.

## Claim

(b) Convention: `removeVertex` contracts through the deleted vertex; deleting i = 7 from the relation 7..10 of 00000030 gives a length-3 relation on 6..9, i.e. 00000300. The pair 90000000 (kA_10/rad^9) and 00000030 (E-190: equivalent, D_10 class, not joined by orbit+mirror moves). The rule of E-114/E-206 (delete a vertex at D >= 3 from every relation support) **does not apply to 90000000**: its single relation covers all ten vertices, so D = 0 everywhere and there is no D >= 3 cell. So the pair cannot test the free-stretch rule. Without the D filter: every deletion of 90000000 gives 80000000 (n = 9, kA_9/rad^8), class `P^(6)_(1,1)` (145 of 1430 rows). 00000030 gives 00000300 (relation 6..9) for i = 1..7, again `P^(6)_(1,1)`, and 0000002 (relation 7..9) for i = 8..10, the hereditary class `P^(0)_(0,8)` (128 rows). Its D >= 3 cells (i = 1..4) land in `P^(6)_(1,1)`, the class of every deletion of 90000000. So (a fact about 00000030's own deletions, 7 of 10, not a transport rule): the images of the two members share a class label for i = 1..7 on 00000030 and every i on 90000000; they split only when the deleted vertex is inside the end of the 3-relation (i = 8, 9, 10, D = 0). This is exact at the n = 9 level; it says nothing about orbit-plus-mirror joining of 80000000 and 00000300 at n = 9 (not tested), nor about any other pair.

(a) Exact `classifyLength(10)`: not run. Seeding by the quipu theorem plus move orbits places 2967 of 4862 rows (61.0 %) in 36 classes in 311 s. The 1895 unplaced rows are searched at depth <= 5 by the decay rule; 6 timed searches took 0.9 to 5.7 s (at depth 5, 4, 3, 3, 3, 3). That bounds the main pass at about 1895 x 6 s = 3 h in the worst case, and much less since each search places the rows it reaches. The step I could not size is `resolveMergeCandidates` (depth-6 searches from every member of each Coxeter-collision group, both directions), and E-032 records the earlier engine at 1010 rows in 9 h. So the cost of the exact run is unknown (the 3 h figure is arithmetic on 6 evenly spaced samples, not a bound); the wrapper below runs it in budgeted pieces. Even the plan step takes 311 s, so a full run exceeds 10 minutes.

## Evidence

Image table (script `experimentalist_pair.py`, labels the exact n = 9 cache; `D` as in E-206):

| source | deleted i | D | n = 9 image | class |
|---|---|---|---|---|
| 90000000 | 1..10 | 0 | 80000000 | P^(6)_(1,1) |
| 00000030 | 1, 2, 3, 4 | 6, 5, 4, 3 | 00000300 | P^(6)_(1,1) |
| 00000030 | 5, 6, 7 | 2, 1, 0 | 00000300 | P^(6)_(1,1) |
| 00000030 | 8, 9, 10 | 0 | 0000002 | P^(0)_(0,8) |

Image classes by threshold: D >= 0: 90000000 one, 00000030 two, union two; D >= 1, 2, 3: 90000000 none, 00000030 one (union one).

Label method and limits. Labels at n = 9 are exact (classifyLength, all 1430 rows, no candidates, cached). The source pair needed no key+orbit label: their equivalence is E-190's Cartan/D_10 argument, a fact about n = 10, not a label. The key+orbit proxy of E-206 is not used here at all. Limit: the n = 9 side being exact does not make the claim general; it is one pair, and the pair is degenerate (a rad^9 LNA has no free stretch). The nearest informative reading is that the "same image" result for i = 1..7 is just 00000030's own deletion behaviour: 7 of 10 deletions give one class, so the agreement with 90000000 is not evidence for a transport rule.

Proposal for `OVERNIGHT.md` (Menu 4), exact labels at n = 10 for E-206's rerun:

```
for k in 1 2 3 4 5 6; do timeout 9.5h .venv/bin/python workshop/rounds/073/experimentalist_plan10.py --run --n 10 --budget-hours 9 --out logs/A_10_exact.csv; rc=$?; [ $rc -eq 2 ] || [ $rc -eq 124 ] || break; done
```

`--run` calls `classifyLength(10, 6, 6, out, resume=True, budgetSeconds=...)` and exits 2 when the budget stopped it, 0 when done; `--budget-hours` is only checked between classes, so one long class step can make `timeout` kill the pass (exit 124): the loop treats 124 like 2 and resumes from the CSV. Worst case 6 x 9.5 h = 57 h; exit 0 or any other code stops the loop (the same contract as `classify.py 10 --resume --budget-hours H`, which already exists; the script only adds the plan). When finished, the exact classes replace the key+orbit labels in `workshop/rounds/070/experimentalist_free10.py` (change `lab` to read the CSV): this removes the 176 unresolved LNAs and the 32-of-42 unfitted caveat of E-206. Resume is by the CSV and its progress sidecar. Output CSV about 0.5 MB: not to be committed.

## Reproduction

```
.venv/bin/python workshop/rounds/073/experimentalist_pair.py                       # 2 s
timeout 10m .venv/bin/python workshop/rounds/073/experimentalist_plan10.py --samples 6 --sample-seconds 60   # about 5.5 min (seed step 311 s)
```

## Prior record

E-190 (the pair, D_10 class), E-114/E-146/E-200/E-206 (free-stretch rule, D >= 3), E-032 and E-008 ("do not schedule `classify.py 10` as a whole-range run", 1010 rows in 9 h on the old engine). New: the pair is outside the rule's domain, and the image classes; the seed-step numbers (61.0 %, 36 classes, 311 s) and search timings at the current engine. Not found in `research/`: a count of how many n = 10 rows the seed step places.

## Code changed

None in the package. New: `workshop/rounds/073/experimentalist_pair.py`, `experimentalist_plan10.py`. No tests (no package code touched).

## Next

Chair: if someone has a night, run the loop above (see the n = 10 time range). Theorist: 80000000 and 00000300 at n = 9 share a class; is that pair joined by orbit+mirror moves (the E-190 phenomenon one length down)? Experimentalist: with exact n = 10 labels, redo E-206 including the 32 unfitted labels and the 176 dropped LNAs.
