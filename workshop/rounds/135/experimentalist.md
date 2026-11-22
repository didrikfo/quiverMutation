# On five lone-3 key classes not used in E-325 (n = 13 h = 1, 4; n = 15 h = 2, 4, 5; head end, key level) no move edge between two K >= 1 rows changes the image key. Every key change on the K=0 / K>=1 boundary is a directed `rules` rewrite K = 0 -> K = 1 onto a row starting 03 (none in reverse), and there are none at n = 13 h = 4 or n = 15 h = 4, 5

author: experimentalist · round: 135 · kind: result (revised after skeptic's review)
thread: S-1 · bears on: H-026, E-325, E-307, E-301
scope: odd n = 13, 15 only; lone-3 Coxeter-key classes n = 13 h = 1, 4 and n = 15 h = 2, 4, 5 (whole class, every K, no cap; 16 439 / 5 023 / 35 745 / 21 045 / 18 328 rows), plus E-325's n = 13 h = 2, 3 and n = 15 h = 3 for the boundary census; head end only; boundary edges counted by direction (the move generator is not symmetric); move families rules / mut / der as in E-325 (library moves assumed derived, not re-checked with the engine); key level; nothing at even n, n >= 17, or for non-lone classes
power: the same edge census with a planted label, the tail-deletion key key(n-1, m[:-1]), finds key-changing edges with both ends K >= 1 in every run: 68 to 4 972 per family and class, against 0 for the head image. A random relabelling finds 1 771 to 89 740. The test can say NO, and does say NO on a label for which no induced-move argument exists.

## Response to referee

Skeptic, round 135, minor revision. All four required items done.

1. **(b) corrected.** The "majority reaches the largest K" half is now stated for `der` bridges only. The two `mut` exceptions are named in (b) and printed in Table B: n = 13 h = 4, 23 rows, {0: 2, 1: 4} with class max K = 5; n = 15 h = 4, 123 rows, {0: 4, 1: 4} with class max K = 7. Read from the existing outputs. I also printed the max-K maps of the other non-largest bridges, which the table had hidden. A third small `mut` bridge has a weak minority: n = 15 h = 5, 25 rows, {0: 6, 1: 2}. It is not an exception, because the majority key reaches the class max K = 6. "Follows from E-301 and E-304" is dropped. Those entries are about whole-class K-tables and do not imply anything bridge by bridge.
2. **Prior record.** It now cites E-325's Power line "NO on all 8 / 16 / 32 key-changing edges". Those are (a)'s counts at n = 13 h = 3, n = 13 h = 2 and n = 15 h = 3. The "not in `research/`" line is narrowed to the shapes, the direction, h = 1, n = 15 h = 2 and the zeros at n = 13 h = 4 and n = 15 h = 4, 5.
3. **Direction.** New script `experimentalist_boundary_dir.py` counts key-changing moves with exactly one end at K = 0, in each direction, under `mut` and `der`. Result: K = 0 -> K >= 1 gives 392 / 16 / 8 / 0 at n = 13 h = 1-4. K >= 1 -> K = 0 gives 0 in every class. Every move is `rules`, the free end has K = 1 and m[4:] is unchanged. Shapes and counts are identical under `mut` and `der`. At n = 15 h = 2, 3 the counts are 64 (3404 -> 0350) and 32 (3405 -> 0360), with 0 in reverse. This reproduces the referee's check. It also shows the `der` half at n = 15, which I had asserted without showing. (a) now writes the shapes with -> and states the asymmetry. The union-find bridge analysis is undirected and unchanged.
4. **Wording.** The "at h >= 4" sentence, the title and *What would refute it* now name n = 13 h = 4 and n = 15 h = 4, 5. The reproduction block has one call per h.

Unchanged: the main claim (0 both-K >= 1 key changes on all six tables) and Table A.

## Claim

On each of the five new classes, under each of `rules`, `mut` and `der`, the number of move edges that change the head-deletion image key key(n-1, m[1:]) with both ends at K >= 1 is **0**. So E-325's claim holds on these orbits: the image key changes only along stretches through K = 0 rows, and restricted components of different image keys meet only through K = 0 "bridge" components.

The new orbits include the central classes n = 13 h = 4 (E-141's 4349-orbit) and n = 15 h = 5 (E-263's class). Neither had its bridges examined before.

Two sharpenings, observed and not proved:
- **(a) Boundary edges.** Every key-changing move with one end at K = 0 and the other at K >= 1 is a `rules` rewrite **directed K = 0 -> K = 1**; the reverse direction (a move from a K >= 1 row to a K = 0 row that changes the key) occurs 0 times, under `mut` and `der`. The K = 1 end starts 03. The head shapes (K = 0 row -> K = 1 row) are:
  - h = 1: four shapes, 4000 -> 0300 (256 edges), 3400 -> 0300 (64), 4030 -> 0302 (64), 3403 -> 0340 (8);
  - h = 2: 3404... -> 0350... only;
  - h = 3: 3405... -> 0360... only.
  In every case the tail m[4:] is unchanged and free. Across 8 classes and n = 13, 15 there are 392, 16, 8, 0 such directed edges at n = 13 h = 1, 2, 3, 4, and 64, 32, 0, 0 at n = 15 h = 2, 3, 4, 5. So at n = 13 h = 4 and n = 15 h = 4, 5 (the only h >= 4 classes checked) every key change lies on an edge between two K = 0 rows.
- **(b) Bridge depth.** In every `mut` and `der` bridge, the minority image keys reach at most K = K0. In every `der` bridge listed, the majority key also reaches the largest K in the class. This second half fails for two small `mut` bridges: n = 13 h = 4, 23 rows (`20270000300`), max K {0: 2, 1: 4} with class max K = 5; n = 15 h = 4, 123 rows (`2000270000300`), {0: 4, 1: 4} with class max K = 7. Under `der` both merge into larger bridges. This is an observation on these bridges only; E-301 and E-304 are statements about whole-class K-tables and do not imply it bridge by bridge.

What this does not claim:
- That the 0 is a theorem. E-325 (ii) says it is expected to be definitional (one-point-extension restriction), so the five new zeros are a code check on new data, not new mathematics.
- Any law for the boundary shapes in (a) beyond the 8 classes.
- Anything about derived classes, as opposed to orbits.

What would refute it: one key-changing edge with both ends K >= 1 in any lone-3 class; for (a), a boundary key change at n = 13 h = 4 or n = 15 h = 4, 5 (already ruled out on these) or, as a test of the pattern, at n = 17 h >= 4.

## Evidence

Table A: key-changing directed edges, grouped by the K of their two ends. The counts are the same for `mut` and `der`; `rules` counts are a subset.

| class | rows | K0 | image: both >= 1 | image: boundary | image: both 0 | tail control: both >= 1 (mut) | random: both >= 1 (mut) |
|---|---|---|---|---|---|---|---|
| n=13 h=1 | 16 439 | 1 | 0 | 392 | 6 320 | 2 952 | 27 102 |
| n=13 h=3 (E-325, re-run) | 5 835 | 3 | 0 | 8 | 1 630 | 783 | 8 425 |
| n=13 h=4 (central) | 5 023 | 4 | 0 | 0 | 1 348 | 642 | 6 725 |
| n=15 h=2 | 35 745 | 2 | 0 | 64 | 10 098 | 4 972 | 77 794 |
| n=15 h=4 | 21 045 | 4 | 0 | 0 | 4 344 | 2 132 | 40 629 |
| n=15 h=5 (central) | 18 328 | 5 | 0 | 0 | 3 212 | 1 566 | 34 078 |

The `rules` both >= 1 count is 0 in all 6 rows. The tail control under `rules` gives 70 to 1 520.

Table B: mixed orbits (`mut`) and their bridges. "Max K per side" gives, for each image key, the largest K that key reaches in the restricted components the bridge touches.

| class | mixed orbit, rows | K0 rows by image key | restricted comps | K = 0 comps | bridges | largest bridge, rows; max K per side |
|---|---|---|---|---|---|---|
| n=13 h=1 | 16 389 | 3141/257/135/173/75 (5 keys) | 33 | 9 | 1 | 10 488; {0: 8, 1-4: 1} |
| n=13 h=4 | 4 349 | 38/32 | 9 | 5 | 3 | 1 791; {0: 5, 1: 4} (others 762 {0: 5, 1: 4}; 23 {0: 2, 1: 4}) |
| n=15 h=2 | 35 505 | 3216/512 | 7 | 3 | 1 | 21 093; {0: 9, 1: 2} |
| n=15 h=4 | 19 963 | 285/128 | 11 | 5 | 3 | 8 751; {0: 7, 1: 4} (others 2 629 {0: 7, 1: 4}; 123 {0: 4, 1: 4}) |
| n=15 h=5 | 14 454 | 64/65 | 10 | 6 | 3 | 4 449; {0: 6, 1: 5} (others 3 574 {0: 6, 1: 5}; 25 {0: 6, 1: 2}) |

Under `der` the bridge count drops to 1, 2, 1, 2, 2. Under `rules` there are no bridges at h = 4 and 5 (0 of 1 804 `rules` K = 0 components at n = 13 h = 4). At h = 1 there are 6 bridges, and at n = 15 h = 2 there are 7. At h = 1, 2 the minority side of a `rules` bridge reaches K0. This answers E-325's open item: at n = 15 h = 2 the minority side reaches K0 = 2, the same as n = 13 h = 2. With h = 3 at K = 1 < K0 (E-325), this still does not separate depth from K0.

The boundary census (`experimentalist_boundary.py`, `mut`, K = 0 side only) and its directed re-run under both `mut` and `der` (`experimentalist_boundary_dir.py`, both directions):
- n = 13: h = 1 gives 392 edges (the K = 0 row starts with 3 in 72 of them and with 4 in 320); h = 2 gives 16; h = 3 gives 8.
- n = 15: h = 2 gives 64; h = 3 gives 32.
- All are `rules` edges K = 0 -> K = 1 with m[4:] unchanged; the K >= 1 -> K = 0 direction gives 0 key changes at every class. Examples: 34050000000 -> 03600000000; 3404000000000 -> 0350000000000.
- Counts and shapes are identical under `mut` and `der` (directed re-run: n = 13 h = 1-4 and n = 15 h = 2, 3; the n = 15 h = 4, 5 zeros come from the `mut` bridge run, where the boundary column is 0 under every family).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_bridge.py 13 1   # 43 s   (also 13 3: 21 s, 13 4: 19 s)
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_bridge.py 15 2   # 4 min 17 s (15 4: 3 min 24 s, 15 5: 3 min 30 s)
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary.py 13 1   # ~15 s
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary.py 13 2   # ~15 s
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary.py 13 3   # ~15 s
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary.py 15 2   # ~3 min
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary.py 15 3   # ~3 min
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary_dir.py 13 H   # one call per H = 1, 2, 3, 4; mut and der, both directions
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary_dir.py 15 H   # one call per H = 2, 3
```

Outputs are in `workshop/rounds/135/experimentalist_bridge_n*_out.txt`, `experimentalist_boundary_n1{3,5}_out.txt` and `experimentalist_boundary_dir_n1{3,5}_out.txt`, all under 4 KB. On this box the key scan is 90% of the n = 15 time.

## Prior record

- E-325 states the claim on n = 13 h = 2, 3 and n = 15 h = 3, and calls (ii) definitional. This submission extends it to five more classes; it is not new mathematics.
- E-141 and E-263 found the mixed orbits of the two central classes. E-301 and E-304 fit K0 on whole-class K-tables; they do not give (b) bridge by bridge.
- E-325's Power line already has the boundary counts at n = 13 h = 3, 2 and n = 15 h = 3: "NO on all 8 / 16 / 32 key-changing edges". These match (a)'s 8, 16, 32. Not in `research/`: the head shapes, the direction, the h = 1 and n = 15 h = 2 counts, and the zeros at n = 13 h = 4, n = 15 h = 4, 5.
- The round-125 Next item, bridge depth at n = 15 h = 2, is answered above.
- Nothing in RETRACTIONS bears on this.

## Code changed

None. New scripts: `workshop/rounds/135/experimentalist_bridge.py`, `experimentalist_boundary.py` and (revision) `experimentalist_boundary_dir.py`. They import `rounds/123/experimentalist_classlevel.py` unchanged.

## Next

- Theorist: prove the restriction property (moves between K >= 1 rows are image moves). That turns Table A's zeros into a lemma. Then explain (a): why the only boundary key changes are these few `rules` shapes onto 03.., and why it needs h <= 3.
- Experimentalist: test whether (a)'s zero at h >= 4 extends to n = 17 h = 4, 6, 7. This needs the faster key or an overnight run, since the n = 17 key scan is about 8x n = 15.
- Experimentalist: test whether the boundary counts double with each extra free tail vertex. 8 -> 32 at h = 3 is a factor of 4 over 2 more vertices, and 16 -> 64 at h = 2 is the same.
