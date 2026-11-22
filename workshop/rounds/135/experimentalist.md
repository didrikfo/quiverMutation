# On five lone-3 key classes not used in E-325 (n = 13 h = 1, 4; n = 15 h = 2, 4, 5; head end, key level) no move edge between two K >= 1 rows changes the image key. Every key change on the K=0 / K>=1 boundary is a `rules` rewrite onto a K = 1 row starting 03, and there are none at h >= 4

author: experimentalist · round: 135 · kind: result
thread: S-1 · bears on: H-026, E-325, E-307, E-301
scope: odd n = 13, 15 only; lone-3 Coxeter-key classes n = 13 h = 1, 4 and n = 15 h = 2, 4, 5 (whole class, every K, no cap; 16 439 / 5 023 / 35 745 / 21 045 / 18 328 rows), plus E-325's n = 13 h = 2, 3 and n = 15 h = 3 for the boundary census; head end only; move families rules / mut / der as in E-325 (library moves assumed derived, not re-checked with the engine); key level; nothing at even n, n >= 17, or for non-lone classes
power: the same edge census with a planted label, the tail-deletion key key(n-1, m[:-1]), finds key-changing edges with both ends K >= 1 in every run: 68 to 4 972 per family and class, against 0 for the head image. A random relabelling finds 1 771 to 89 740. The test can say NO, and does say NO on a label for which no induced-move argument exists.

## Claim

On each of the five new classes, under each of `rules`, `mut` and `der`, the number of move edges that change the head-deletion image key key(n-1, m[1:]) with both ends at K >= 1 is **0**. So E-325's claim holds on these orbits: the image key changes only along stretches through K = 0 rows, and restricted components of different image keys meet only through K = 0 "bridge" components.

The new orbits include the central classes n = 13 h = 4 (E-141's 4349-orbit) and n = 15 h = 5 (E-263's class). Neither had its bridges examined before.

Two sharpenings, observed and not proved:
- **(a) Boundary edges.** Every key-changing edge with one end at K = 0 and the other at K >= 1 is a `rules` rewrite whose free end has K = 1. Its free end starts 03. The head shapes of the K = 0 end and the free end are:
  - h = 1: four shapes, 4000 <-> 0300 (256 edges), 3400 <-> 0300 (64), 4030 <-> 0302 (64), 3403 <-> 0340 (8);
  - h = 2: 3404... <-> 0350... only;
  - h = 3: 3405... <-> 0360... only.
  In every case the tail is free. Across 8 classes and n = 13, 15 there are 392, 16, 8, 0 such edges at n = 13 h = 1, 2, 3, 4, and 64, 32, 0, 0 at n = 15 h = 2, 3, 4, 5. At h >= 4 every key change lies on an edge between two K = 0 rows.
- **(b) Bridge depth.** In every `mut` and `der` bridge, the minority image keys reach at most K = K0, and the majority key reaches the largest K in the class. This is H-026 read on the bridges, so it follows from E-301 and E-304 and is not independent evidence.

What this does not claim:
- That the 0 is a theorem. E-325 (ii) says it is expected to be definitional (one-point-extension restriction), so the five new zeros are a code check on new data, not new mathematics.
- Any law for the boundary shapes in (a) beyond the 8 classes.
- Anything about derived classes, as opposed to orbits.

What would refute it: one key-changing edge with both ends K >= 1 in any lone-3 class, or a boundary key change at h >= 4.

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
| n=13 h=4 | 4 349 | 38/32 | 9 | 5 | 3 | 1 791; {0: 5, 1: 4} (others 762, 23) |
| n=15 h=2 | 35 505 | 3216/512 | 7 | 3 | 1 | 21 093; {0: 9, 1: 2} |
| n=15 h=4 | 19 963 | 285/128 | 11 | 5 | 3 | 8 751; {0: 7, 1: 4} (others 2 629, 123) |
| n=15 h=5 | 14 454 | 64/65 | 10 | 6 | 3 | 4 449; {0: 6, 1: 5} (others 3 574, 25) |

Under `der` the bridge count drops to 1, 2, 1, 2, 2. Under `rules` there are no bridges at h = 4 and 5 (0 of 1 804 `rules` K = 0 components at n = 13 h = 4). At h = 1 there are 6 bridges, and at n = 15 h = 2 there are 7. At h = 1, 2 the minority side of a `rules` bridge reaches K0. This answers E-325's open item: at n = 15 h = 2 the minority side reaches K0 = 2, the same as n = 13 h = 2. With h = 3 at K = 1 < K0 (E-325), this still does not separate depth from K0.

The boundary census uses `mut` and is identical under `der`:
- n = 13: h = 1 gives 392 edges (the K = 0 row starts with 3 in 72 of them and with 4 in 320); h = 2 gives 16; h = 3 gives 8.
- n = 15: h = 2 gives 64; h = 3 gives 32.
- All are `rules` edges with the free end at K = 1. Examples: 34050000000 -> 03600000000; 3404000000000 -> 0350000000000.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_bridge.py 13 1   # 43 s   (also 13 3: 21 s, 13 4: 19 s)
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_bridge.py 15 2   # 4 min 17 s (15 4: 3 min 24 s, 15 5: 3 min 30 s)
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary.py 13 3 # ~15 s each for h = 1, 2, 3
timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary.py 15 2 # ~3 min each for h = 2, 3
```

Outputs are in `workshop/rounds/135/experimentalist_bridge_n*_out.txt` and `experimentalist_boundary_n1{3,5}_out.txt`, all under 4 KB. On this box the key scan is 90% of the n = 15 time.

## Prior record

- E-325 states the claim on n = 13 h = 2, 3 and n = 15 h = 3, and calls (ii) definitional. This submission extends it to five more classes; it is not new mathematics.
- E-141 and E-263 found the mixed orbits of the two central classes. E-301 and E-304 fit K0, which is the source of (b).
- The boundary-shape census (a) and "no boundary key change at h >= 4" are not in `research/`. I grepped for "boundary", "bridge" and "K = 0".
- The round-125 Next item, bridge depth at n = 15 h = 2, is answered above.
- Nothing in RETRACTIONS bears on this.

## Code changed

None. New scripts: `workshop/rounds/135/experimentalist_bridge.py` and `experimentalist_boundary.py`. They import `rounds/123/experimentalist_classlevel.py` unchanged.

## Next

- Theorist: prove the restriction property (moves between K >= 1 rows are image moves). That turns Table A's zeros into a lemma. Then explain (a): why the only boundary key changes are these few `rules` shapes onto 03.., and why it needs h <= 3.
- Experimentalist: test whether (a)'s zero for h >= 4 holds at n = 17 h = 4, 6, 7. This needs the faster key or an overnight run, since the n = 17 key scan is about 8x n = 15.
- Experimentalist: test whether the boundary counts double with each extra free tail vertex. 8 -> 32 at h = 3 is a factor of 4 over 2 more vertices, and 16 -> 64 at h = 2 is the same.
