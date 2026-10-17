# At n = 10 -> 9 the free-stretch rule (distance >= 3 from every relation) gives one image class for 10 of 10 key+orbit labels that have such a cell; 32 of 42 labels have none, and D >= 2 fails for 2 of 19 labels, both inside a single exact class (extends E-200 from n = 9 to n = 10; new only in the D-cell counts and the unfitted list)

author: experimentalist · round: 070 · kind: rediscovery
thread: S-1 · bears on: E-114, E-146, E-194, E-200, H-020
scope: n = 10 -> 9 only; 4 686 of 4 862 LNAs (176 cospectral-unresolved dropped, as in E-114); source labels are Coxeter key + free/edge/double orbit + mirror (not `classifyLength`: a full n = 10 run takes hours; single-class status checked only for the two D = 2 failing labels, see Response); image labels are exact `classifyLength` classes at n = 9 (round-065 cache, 20 classes); derived classes only; orbit, verdict, head/tail not tested; one deleted vertex per cell; no search over other rules.

## Response to referee

1. Cell count: corrected. Three labels rest on 6 or fewer D >= 3 cells (104 rows: 2; 78 rows: 6; 54 rows: 6), so "10 of 10" is thin in the same way as E-200's 6 of 6 (4 and 2 cells). The 382 cells are not 382 independent checks; they test "one image per label", nothing more.
2. Pooling: the caveat was pointed the wrong way and is replaced (see Claim). Pooling cannot hide a D >= 3 failure; it can manufacture D = 1 and D = 2 failures. Tested directly with `experimentalist_pool.py` (about 2 min, output `experimentalist_pool.txt`): the 262-row label is one free+edge+double orbit joined with its mirror, hence one exact class (orbit moves are known equivalences). The 680-row label is three such components (368, 242, 70 rows); `meetingPoints` (depth 2 to 4, with duals) between their representatives found no link, which proves nothing, so I do not claim the label is one class. But the 368-row component alone, a single class, already has both images at D >= 2 (P^(5)_(1,2) 32 cells, P^(4)_(1,3) 100). So both D = 2 failures are real failures inside one exact class, not artefacts of pooling. "Not by relation length" is kept: inside the 368-row component, next to a 3-relation the image is P^(5)_(1,2) in 16 cells and P^(4)_(1,3) in 11 (same split on the left and right); the 262-row label is a single class. The earlier 16 vs 15 count pooled the three components.
3. E-200 is now cited as the n = 9 -> 8 analogue (free-stretch d >= 3: 6 of 6 classes; d >= 2 fails in P^(4)_(1,3)). This entry extends it to n = 10 and adds the D-cell table and the unfitted list.
4. 5 326 / 5 326 dropped from the claim; it is a consequence of 10 / 10 (one image per label forces all same-label pairs to agree) and adds no evidence.
5. Title and text say "labels" for key+orbit labels.

## Claim

At n = 10 -> 9, deleting a vertex at distance D >= 3 from the support of every relation (382 (row, i) cells, all on a free end: 0 gap cells) sends each of the 10 key+orbit labels that have such a cell to a single n = 9 class; for D >= 4 it is 5 of 5 (126 cells). This is the n = 10 analogue of E-200's n = 9 -> 8 result. Three of the ten labels rest on 6 or fewer cells (2, 6, 6). The rule is silent for 32 of the 42 labels (no vertex at D >= 3), among them the large ones (312, 166, 164, 123, 111 rows ...). D >= 2 fails for 2 of 19 labels (680 rows and 262 rows, 2 images each), and both failures occur inside a single exact class (the 262-row label is one class; the 680-row label contains a 368-row single-class component that fails by itself). D >= 1 fails for 13 of 34 labels; those failures may be inflated by pooling. Pooling cannot hide a D >= 3 failure (a class with two images makes every label containing it have two images), so a full `classifyLength(10)` labelling can only split labels and cannot turn the 10 of 10 into a failure; it could only make some labels' cell counts smaller. Refuted by: a label (or class) with a D >= 3 cell and two image classes (found none).

## Evidence

| D >= | labels with a cell | one image | not | cells | labels with no cell |
|---|---|---|---|---|---|
| 1 | 34 | 21 | 13 | 4 124 (184 interior-gap) | 8 |
| 2 | 19 | 17 | 2 | 1 204 (8 gap) | 23 |
| 3 | 10 | 10 | 0 | 382 (0 gap) | 32 |
| 4 | 5 | 5 | 0 | 126 | 37 |

- D >= 3 images, source label (key tail) -> image: (..-1,-1,-1,-1,0,1,1) 1125 rows -> P^(5)_(1,2) (88 cells); (..-2,-2,-2,-1,0,1,1) 680 -> P^(4)_(1,3) (44); 323 rows -> P^(6)_(1,1) (64); 313 -> P^(4)_(1,3) (18); 262 -> P^(3)_(2,3) (12); 256 (hereditary) -> P^(0)_(0,8) (124); 165 -> P^(4)_(2,2) (18); 104 -> P^(1,1)_(1,3,1) (2); 78 -> P^(1,3)_(1,1,1) (6); 54 -> P^(2,3)_(1,0,1) (6). Three labels rest on 6 or fewer cells (2, 6, 6). Full list: `experimentalist_free10.txt`.
- The D >= 2 failures are the same two labels as in E-114's K0 = 2 count for n = 10 (2 of 12 there). In both, the image splits between two classes by the vertex, and not by the length of the nearest relation: the 262-row label (one exact class) next to a 4-relation gives P^(4)_(2,2) 8 and P^(3)_(2,3) 9 (`experimentalist_free10_d2.txt`); in the 368-row single-class component of the 680-row label, next to a 3-relation P^(5)_(1,2) 16 and P^(4)_(1,3) 11 (`experimentalist_pool.txt`; the pooled-label count 16 vs 15 in `_d2.txt` mixes three components). This answers the "does the first relation's length decide" question of my round-067 notebook: no, at n = 10.
- No D >= 3 cell lies in a gap between relations at n = 10 (a gap of 5 plus two supports of 2 does not fit with a core); the 8 interior cells at D = 1 and D = 2 are too few to test gap deletion alone.
- Counts per label reproduce E-146's table for n = 10 (K0 = 3: 12 of 12 there, 10 of 10 here: E-146 counts classes with an end of length >= 3, mine counts vertices at distance >= 3 from supports, which excludes the end vertices adjacent to a relation; the two differ by definition, not by data).

## Reproduction

```
.venv/bin/python workshop/rounds/070/experimentalist_free10.py > workshop/rounds/070/experimentalist_free10.txt     # 10 s
.venv/bin/python workshop/rounds/070/experimentalist_free10_d2.py > workshop/rounds/070/experimentalist_free10_d2.txt # 12 s
timeout 10m .venv/bin/python workshop/rounds/070/experimentalist_pool.py > workshop/rounds/070/experimentalist_pool.txt   # about 2 min
```
Needs `workshop/rounds/027/maverick_classes.py` and `workshop/rounds/065/experimentalist_classes_n9.csv`. `--plan` was not needed: the sizing is 4 862 rows x 10 deletions, no search.

## Prior record

E-114 (K >= 3 rule at n = 8, 9, 10: 3, 5, 10 classes), E-146 (corrected labels: n = 10 12/12 at K0 = 3, 2 failures at K0 = 2), E-200 (same at n = 9 -> 8: free-stretch d >= 3 gives 6 of 6, d >= 2 fails in P^(4)_(1,3); this entry is its n = 10 analogue), E-173 (n = 15 sized). This result is a rediscovery at n = 10; new only: the distance-D cell table, the explicit unfitted list (32 of 42), the single-class check of the two D = 2 failures, and the no-by-relation-length answer for the D = 2 splits. It is not a new test of a law, since the n = 10 labels are the ones E-114 used.

## What does not fit in 10 minutes (proposal for `OVERNIGHT.md`)

1. `classifyLength(10)` with `resume=True`, exact classes for all 4 862 + 176 LNAs (hours, per the docstring), then rerun `experimentalist_free10.py` with the exact labels to remove the key+orbit caveat and place the 176 unresolved rows. 
2. The unfitted n = 10 labels (32) need a different rule: n = 11 -> 10 (16 796 rows, 442 unresolved in E-114's labels) would add D >= 3 cells for more classes; keys only, about a minute, but E-114 notes n = 11 label gaps. Not run, I did not size it with `--plan`.
3. E-173's n = 15 class-level scan (2.67M rows, ~33 min on 4 cores, 24 shards of ~6 min): can be run as 24 foreground-sized jobs; not run, out of scope for a single 10-minute command and unchanged in expectation.

## Code changed

None (two new scripts only). No tests needed.

## Next

Experimentalist: item 1 above once the human opens the overnight slot. Theorist: why D >= 3 is the threshold (F-028 free move?) and what splits the two D = 2 failures; a rule for the 32 unfitted labels. S-1 items 2 (head/tail, orbit) remain open at n >= 12 (round-067 notebook item 1).
