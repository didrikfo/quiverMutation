# At n = 9 -> 8, a full-20-class recount of E-114's fixed-vertex deletion rates (same numbers, 37.5-48.9%); no fixed i makes every class pure, and some (i, i') agrees for every pair (trivial)

author: experimentalist · round: 065 · kind: negative
thread: S-1 · bears on: F-011, F-012, E-146, E-173
scope: refines E-114 (not new): n = 9 (all 1430 rows, 20 classes, F-011) to n = 8 (429 rows, 11 classes); fixed vertex i or fixed pair (i, j) only, not a rule depending on the row (E-114's free-end K >= 3 rule is one); classes by `classifyLength` (resolves the 16 rows E-114 dropped); no sampling.

## Response to referee

1. Done. E-114 and synthesis 001 section 4 cited in Prior record; "no entry counts" and "this is new" deleted; title and Scope now say
   refinement. Differences stated: class source (key + orbit vs `classifyLength`, 16 unresolved rows now resolved), unordered vs
   ordered pairs (147 681 vs 295 298), and what is new.
2. Done. 0.4887 / 0.4285 agree with E-114's 0.489 / 0.429 to the digits shown; range and baseline lie inside its ranges; no differing number.
3. Done. 9 x 9 matrix in the Evidence section; output in `experimentalist_delcompat_n9.txt`; the /tmp CSVs are over 200 KB so the
   derived (relation, class) maps are `experimentalist_classes_n8.csv`, `_n9.csv`; regeneration in Reproduction.
4. Done. "Refuted by" now says fixed i or (i, j); E-114's free-end K >= 3 rule is named as a row-dependent rule not excluded here, and why
   it does not contradict "no fixed i".
Not done: independent cold re-classification at n = 9 (over 10 min per call; the referee also took the classes as given). Row-by-row
comparison with E-114's labels: its row-level output is not in the repository, so I state only the totals.

## Claim

S-1 item 1, counted exhaustively at n = 9. There are 147 681 same-class pairs of LNAs at n = 9. For a fixed
deleted vertex i applied to both members (`piecewiseHereditary.removeVertex`), the fraction whose images lie in
one n = 8 class is: i = 1,2,8,9: 0.4887; i = 3,7: 0.4471; i = 5: 0.4285; i = 4,6: 0.3750 (baseline for two
random n = 8 rows: 0.176). No i is compatible: for no i does every n = 9 class land in a single n = 8
class (class sizes 377 / 300 / 145 map to 5 / 6 / 2 classes at every i, only the 128-row almost-hereditary
class `P^(0)_(0,8)` and the three small classes map to 1). Allowing different vertices i, i' on the two
members, the best off-diagonal cell is 0.492 (i, i' in {1,2} x {8,9}), still under one half. What it does
NOT claim: nothing about n >= 10, about non-uniform rules (a rule depending on the relation pattern), or about
the Coxeter key (E-146/E-173 are key-level). Every removeVertex image at n = 9 was admissible (0 undefined).
Refuted by: a fixed i (or fixed pair i, j) giving image-class agreement on all pairs. A row-dependent rule i(Lambda)
is not excluded: E-114's free-end rule (delete a vertex of a free end of K >= 3 vertices) is one, and gives an image class
depending only on the source class for the resolved classes that have such an end. It does not contradict "no fixed i",
since the deleted position there changes with the row (the end where the relations are absent).

Existence, by contrast, is trivial: for every one of the 147 681 pairs there is some (i, i') with
delta_i x ~ delta_i' y (the image-class sets always meet), and 94.8% (139 945) already with i = i'. This is
the weak statement and should not be read as compatibility; it says only that the 11 target classes are
few and each n = 9 class reaches a common one.

## Evidence

Fractions of same-class pairs preserved, same i on both members (defined pairs = 147 681 for every i):

| i | 1,2,8,9 | 3,7 | 4,6 | 5 |
|---|---|---|---|---|
| same-class pairs kept | 72 177 | 66 025 | 55 382 | 63 275 |
| fraction | 0.4887 | 0.4471 | 0.3750 | 0.4285 |

Matrix (x deleted at i, y at j, rows i, columns j; same-class unordered-pair fractions):

```
MATRIX fraction (x deleted at i, y at j), rows i, cols j
1 0.489 0.489 0.461 0.411 0.422 0.411 0.463 0.492 0.492
2 0.489 0.489 0.461 0.411 0.422 0.411 0.463 0.492 0.492
3 0.461 0.461 0.447 0.395 0.403 0.397 0.450 0.463 0.463
4 0.411 0.411 0.395 0.375 0.387 0.378 0.397 0.411 0.411
5 0.422 0.422 0.403 0.387 0.428 0.387 0.403 0.422 0.422
6 0.411 0.411 0.397 0.378 0.387 0.375 0.395 0.411 0.411
7 0.463 0.463 0.450 0.397 0.403 0.395 0.447 0.461 0.461
8 0.492 0.492 0.463 0.411 0.422 0.411 0.461 0.489 0.489
9 0.492 0.492 0.463 0.411 0.422 0.411 0.461 0.489 0.489
```

Symmetry i <-> 10 - i is exact (duality) and 1, 2 give identical counts. Full 9 x 9 matrix of (i, j) is below (also in `experimentalist_delcompat_n9.txt`) (range 0.375 to 0.492; larger when the two vertices are both near an end). Per-class image-class
counts (columns i = 1..9), largest classes: `P^(5)_(1,2)` (377 rows) 5 at every i; `P^(4)_(1,3)` (300) 6 at
every i; `P^(6)_(1,1)` (145) 2; `P^(1,3)_(1,1,1)` (44) 4 4 5 5 2 5 5 4 4; `P^(2,3)_(1,0,1)` (30) 2 2 3 4 2 4 3 2 2.
Vertex 5 (centre) is the only i that collapses small classes to 1 image class (`P^(1,2)_(1,1,2)`,
`P^(2,2)_(1,0,2)`, `P^(1,1,2)_(1,0,0,1)`, `P^(1,1,1)_(1,0,1,1)`), but it fails on the large ones
(`P^(5)_(1,2)` still 5).
Table sizes: n = 9 recomputed to 20 classes (matches F-011), n = 8 to 11; first run stopped at 21 classes
because the 10-minute cap hit before the resolve step (3 rows `3450400` unmerged); resumed to 20. A capped
classifyLength is not a verdict, as the persona warns.
Item 2 (predictable change) was not attempted.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/065/experimentalist_delcompat.py 9
```
Classes are computed on first call (n = 8 about 40 s; n = 9 about 15 min in total over 3 resumed calls
because of the cap, 181 s for the last) and cached in `/tmp/s1del_A_<n>.csv`. For a cold start run
`classifyLength(9, ..., resume=True)` on a cold cache, repeatedly under `timeout 10m`
(the script now passes resume=True, so re-run it until it finishes); analysis alone takes 3 s. Output committed: `experimentalist_delcompat_n9.txt` (3 KB: matrix, per-class image counts, existence table). The two /tmp CSVs
(77 KB, 278 KB) are not copied; their (relation string, class) columns are in `experimentalist_classes_n8.csv`, `_n9.csv`
(`experimentalist_dumpclasses.py`, 80 KB total). Regeneration: delete /tmp/s1del_A_*.csv and re-run the command above
until it prints (resume=True continues the cap-interrupted classification).

## Prior record

E-114 (round 027, maverick; synthesis 001 section 4) already reports this count at n = 9: rate 0.37-0.49 for every (i, j) against
0.17-0.19, 0.429 at i = j = middle, 0.489 at (1, 1), "no choice of position", and the free-end K >= 3 rule. Agreement with E-114:
0.4285 and 0.4887 equal its 0.429 and 0.489 to the digits shown, our range 0.375-0.492 sits inside its 0.37-0.49, baseline
0.176 inside 0.17-0.19; no number differs. What differs: (a) class source: E-114 uses the Coxeter key split by the orbit of free, edge,
double moves and the mirror, dropping 16 unresolved LNAs at n = 9; here `classifyLength` gives all 1430 rows in the 20 F-011
classes (the 16 resolved), (b) pairs: E-114 counts 295 298 ordered pairs, here 147 681 unordered (x 2 = 295 362, so the 16 resolved
rows account for 64 ordered pairs; E-114 prints 3 digits), (c) new here: per-class image-class
counts, the best off-diagonal 0.492, the vertex-5 collapse of small classes, existence 100% (94.8% with i = i'). These refine E-114;
I have not compared the two label sets row by row. E-146 / E-173 are key-level; F-012 uses deletion for non-piecewise-heredity
certificates. The earlier text "no entry counts this" was wrong (I grepped S-1 and removeVertex, not the E-114 title).

## Code changed

New script `workshop/rounds/065/experimentalist_delcompat.py`. No library change, no tests run.

## Next

Maverick/theorist: restrict to pairs whose deletion is in a free stretch (item 1's "far from every core")
and recount; my numbers are over all i, with no free-stretch filter. Experimentalist: the same count at
n = 10 -> 9 needs the n = 10 classes (hours; OVERNIGHT proposal), and item 2 (which invariant changes) by
tabulating image class vs source class per i at n = 9 (the 20 x 11 table is a 3 s extension).
