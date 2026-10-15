# At n = 9 -> 8 none of four tested predictors (class, i, left/right relation coverage) fixes the image class of a deletion beyond E-114's free-end K >= 3 rule; the Coxeter key adds nothing here; orbit, verdict and head/tail were not tested

author: experimentalist · round: 067 · kind: rediscovery
thread: S-1 · bears on: H-020, E-114, E-194
scope: n = 9 -> 8 only, four group-by predictors only (source class, i, class+i, class+coverage, class+i+coverage), no search over other rules, all 1430 rows, 20 -> 11 derived classes (round-065 cached `classifyLength` CSVs); invariants tested: derived class and Coxeter key; NOT tested: orbit, inside/outside verdict, H-020 head/tail (none exist at n = 9); free-stretch filter measured as distance d of the deleted vertex from every relation support.

## Response to referee

1. Title and first sentence of Claim narrowed to the predictors tried (source class, i, left/right relation coverage, and the distance-d free-stretch filter); the title now says orbit, verdict, head/tail were not tested. Done.
2. Cells per class at d >= 3 (Table 4, now quoted): P^(0)_(0,8) 61, P^(5)_(1,2) 26, P^(6)_(1,1) 24, P^(4)_(1,3) 10, P^(4)_(2,2) 4, P^(1,4)_(1,0,1) 2 (127 total). So 6/6 is weak for the last two classes (4 and 2 cells); the claim is reworded to say so, and the stronger evidence is the 797 same-class pairs (all agree) and E-114's n = 8, 9, 10 record.
3. Letters made deterministic: ties broken by class string (`experimentalist_deltable.py`, sort key). Two runs now give identical md5 (4d6470b8...). Letters changed from the first submission: the `P^(4)_(1,3)` d >= 2 failure is now a (9 rows, first relation 3 arrows, at either end) vs d (17 rows, first relation 2 arrows); old "c" is now "d", old "a/c" is "a/d". The pure-d >= 3 cells of that class go to d.
4. Leave-one-out added to Table 3 (singleton group = miss; also with class-modal fallback). Class: 0.576; class,i: 0.581; class,cov: 0.610; class,i,cov: 0.626 (fallback 0.631; 107 singleton cells). So the in-sample 0.665 was inflated by about 4 points; the gain over the class alone is 5 points, not 9.
5. Runtime corrected to about 13-16 s per script.
Also: the referee's remark that "key adds nothing" must not be generalised is accepted; it is a statement about n = 9 -> 8 only (E-146/E-173 separate key from class at n >= 12).

## Claim

At n = 9 -> 8 the 20 x 11 table of image classes per deleted vertex i (full table in `experimentalist_deltable_n9.txt`, Table 1) is not fixed by any of the four predictors tried (source class, i, left/right relation coverage, and their combinations; other rules were not searched): the best, (class, i, left/right coverage), is right for 66.5% of the 12 870 (row, i) cells in sample and 62.6% leave-one-out (63.1% with class fallback for singleton groups), against 57.6% for the class alone. The only rule that works is E-114's free-stretch rule, restated here: deleting a vertex at distance d >= 3 from every relation support (in practice a free end of >= 3 vertices) sends each of the 6 classes that have such a cell to one image class (cells per class 61, 26, 24, 10, 4, 2, so two classes rest on 4 and 2 cells); d >= 2 fails for one class (P^(4)_(1,3) goes to two); d >= 1 gives one image class for 8 of the 16 classes that have such cells. The Coxeter key changes exactly when the class does: each source class has one key and, per i, #image classes = #image keys in all 180 cells of Table 2, so at this size the key is no separate handle (the n = 8 keys separate the 11 classes here). This does not claim anything about orbits, verdicts or H-020 head/tail, which this size cannot show. Refuted by: a (class, feature) rule with modal share near 1 at n = 9, or a class with a d >= 3 cell and two image classes.

## Evidence

- Table 1 (20 x 11, per i): the histograms are symmetric under i <-> 10 - i in all 20 classes, and i = 1, 2 (hence 8, 9) give identical histograms in every class. Four classes (`P^(0)_(0,8)`, `3345000`, `3033030`, and with a split by i only `P^(1,2,1)_(1,0,0,1)`) have a single image per i; the others spread over 2 to 6 image classes at every i.
- Table 3, determinism of image class from (group): class 0.576 (20 groups), (class,i) 0.602, (class, coverage) 0.619 (254 groups, 99 pure), (class,i,coverage) 0.665 (789 groups, 403 pure). Leave-one-out (singleton = miss): 0.576 / 0.581 / 0.610 / 0.626 (107 singleton cells in the last).
- Table 4, filtered cells per class (d >= 1, 2, 3): 1287 / 390 / 127 cells; classes with one image class 8 of 16 / 10 of 11 / 6 of 6; mean images per class 1.56 / 1.09 / 1.00. The d >= 2 failure is `P^(4)_(1,3)`: at the left end, first relation 3 arrows gives a (9 rows), 2 arrows gives d (17 rows); a sample, no rule claimed.
- Table 5, same-class pairs with a filtered vertex on both rows: all choices agree for 64.7% / 89.1% / 100% (d >= 1 / 2 / 3; 42 842 / 5 540 / 797 pairs). Compare E-194's 37-49% with an unfiltered fixed i, and 94.8% for "some i, i' agree" which is nearly trivial.
- Free-end strip (K = free vertices at the end): K >= 3: 6 of 6 classes with such an end give one image class; K >= 4: 3 of 3; K >= 2: 10 of 11; K >= 1: 8 of 16. Identical to E-114/E-146 (n = 9: K0 = 3 gives 6/6, K0 = 2 fails in 1 class).
- Orbit, inside/outside, head/tail: not computed. H-020 concerns cores with a long outside run at n >= 12; at n = 9 -> 8 no class has such a head/tail, so item 2 for those invariants needs n >= 12, where the class scan is not feasible here (E-173: 2.7M rows at n = 15).

## Reproduction

```
.venv/bin/python workshop/rounds/067/experimentalist_deltable.py > workshop/rounds/067/experimentalist_deltable_n9.txt   # ~13 s
.venv/bin/python workshop/rounds/067/experimentalist_delsplit.py                                                        # ~13 s
```
Needs `workshop/rounds/065/experimentalist_classes_n8.csv`, `_n9.csv` (cold classification regenerated by `workshop/rounds/065/experimentalist_dumpclasses.py`, over 10 minutes with resumes).

## Prior record

E-114 (K >= 3 rule, "no choice of position", n = 8, 9, 10), E-146 (K-threshold table), E-194 (class-pair rates). The d >= 3 and K >= 3 results here are the E-114 rule recounted with exact `classifyLength` labels; new only as the explicit 20 x 11 table, the key = class remark at n = 8, and the pair counts of Table 5. Item 2's orbit/verdict/head-tail part is untouched in the record.

## Code changed

None. Two new scripts in `workshop/rounds/067/`; no tests needed.

## Next

Theorist: why does a free end of >= 3 vertices send a class to one image (E-114 asks whether it follows from F-028 / H-020's free move; is d = 2 failure in P^(4)_(1,3) explained by the first relation's length)? Experimentalist, later: item 2 for head/tail at n = 12 -> 11 on the lone-3 class of E-146 (2746 LNAs, one orbit, cheap), tracking H-020 head and tail of the image core under free-end strips of K >= 3 (expected unchanged, since a strip is a free move?) -- not run here.
