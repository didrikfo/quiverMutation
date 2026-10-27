# At n = 9, the one key-kept perI J != 0 step in the 2400230 ball (parent depth 8) gives no bridged total <= 13 between 2400230 and 3304000 among the explored balls (C radius 4, 3304000 radius 5, 2400230 radius 6) (null, partial coverage)

author: experimentalist · round: 091 · kind: negative
thread: T10 · bears on: F-010, F-047, E-226, E-238, E-244, E-250
scope: n = 9; E-212 `moves` (gate, J = 0 `tiltingPlus`, key kept) balls: 2400230 to depth 8 (20 210 nodes, all scanned as parents), 3304000 to depth 5 (4 046); J != 0 = the `perI` flag at the gate, key-kept = child key equals the class key; the 3304000 side contributes only the already recorded E-239 / E-244 steps, not rescanned, and E-244's index 8427 was not run; the claim is only for the explored radii (C 4, 3304000 5, 2400230 6), no total <= 13; generation by T assumed.

## Response to referee

Re-run this sitting (shell worked on the 3rd try): the three meets and a new per-depth scan; `/tmp/e91` pickles were still present. The 49 s C-ball and the builds were not rebuilt (pickles reused), so item 6 is only partly met: the meets and scan are reproduced from the saved balls, ball construction is not.
1. Key: `experimentalist_jscan2_out.txt` lines 1-3 print `key(c) == base` for 2400230, 3060000, 3304000: all True (`base` = Coxeter key of `KEYROOT`, default 3060000). So "kept" means the common key of the three classes.
2. Saved: `experimentalist_meets_out.txt` (three meets) and `experimentalist_jscan2_out.txt` (per depth). Correction: perI-flagged parents by depth 0..8 are 0,0,0,0,0,0,3,47,348 (398 parents, 655 steps), so none at depth <= 5 holds, but flags start at depth 6, not 8; kept only at depth 8.
3. The 7 in "1 + 5 of 7": E-239's 1 child plus E-244's 6, disjoint (E-244 is depth-6 parents only); 5 of those 7 join 3304000 by J = 0.
4. Yes: meet C(4) vs 2400230(6) (total 9 = 3 + 6, Hom-verified, iso) is the positive control of the meet tool; reproduced. Re-run results: C-R 0 tests none; C-A 1 test total 9; R-A 0 tests none.
5. Title, claim and Scope narrowed to the explored radii; "gap not closed" now "no bridged total <= 13 among the explored balls".
Also conceded (referee True? 4): the structural "no bridge" holds for any J != 0 step whose child is in its own J = 0 component, independent of the scan; only 8427 and unscanned exits are informative.


## Claim

The calibration asked for ("add the key-kept J != 0 edges of E-234, E-239, E-244 and see whether 2400230 meets 3304000 below 6 + 8") has a structural part and a measured part.

Structural: of the recorded J != 0 children, those of E-234 (3060000 ball) belong to the other F-010 class (3060000 !~ 3304000, E-238, F-010) and cannot bridge to 3304000 at all; those of E-239 and E-244 that join 3304000 by J = 0 (1 + 5 of 7) already lie in 3304000's J = 0 component, so an edge into them only shortens paths inside that component and cannot reach 2400230 unless the component already contains it. The only recorded candidates that could bridge are E-244 index 8427 (not run) and J != 0 steps leaving the 2400230 ball, which no earlier round scanned.

Measured: I scanned all 20 210 nodes (depth <= 8) of the 2400230 J = 0 ball. 655 perI-flagged candidate steps in 398 parents (none at depth <= 5; per depth 0,0,0,0,0,0,3,47,348), 1 is key-kept (parent depth 8, path F1 F1 F2 F3 R9 R7 R9 F3, step F at vertex 3), i.e. 1 in 20 210 parents. Its child C joins the 2400230 ball (iso meet, J = 0, Hom-verified, total 9 = 3 + 6), so it does not leave the 2400230 component, and its radius-4 ball (2 692) has no meet with the radius-5 3304000 ball (4 046): C + 4 + 5, no meet. The 2400230 radius-6 ball (2 691) also has no meet with the 3304000 radius-5 ball: 6 + 5 = 11, no meet (E-226 had 6 + 8 and 5 + 5). So no bridged total at or below 13 exists in the explored sets (a bridge through C needs C-R distance <= 4+5, none); this refutes nothing and confirms nothing about the J = 0 premise. The claim would be refuted by a J != 0 child (from the 2400230 ball, the 3304000 ball or 8427) whose J = 0 ball meets the other root's ball.

Not claimed: that J != 0 steps cannot close the gap in general. The one step found sits at parent depth 8, so any bridge through it costs at least 9 on the 2400230 side, above the 14 of E-226 once the 3304000 side's depth 5 is added; the test as sized could only expose a meet with C within distance 4 of 3304000, and there is none.

## Evidence

| set | size | note |
|---|---|---|
| 2400230 ball, depth 0..8 | 1, 9, 40, 130, 369, 1 001, 2 691, 7 261, 20 210 | build 36+61+174 s etc., about 5 min in two commands |
| parents scanned | 20 210 | 655 perI-flagged key-moved, 1 key-kept (depth 8) |
| C ball radius 4 | 1, 14, 100, 543, 2 692 | 49 s |
| 3304000 ball radius 5 | 4 046 | matches E-226 / E-244 |
| meet C(4) vs 3304000(5) | 0 WL tests | no meet |
| meet C(4) vs 2400230(6) | 1 test | total 9 (3 + 6), Hom-verified, iso |
| meet 3304000(5) vs 2400230(6) | 0 WL tests | no meet (5 + 6) |

The scan rate (1 in 20 210) is a single event and no rate is claimed (E-244 comparison: 6 of 11 355 at 3304000 depth 6, events not independent). That there are no flagged candidates at depth <= 5 in the 2400230 ball differs from 3304000 (354 flagged of 2 958 parents at depth 5, E-239): the 2400230 ball is smaller by a factor 4 at depth 5 (1 001 vs 4 046) and has few perI-capable vertices.

## Reproduction

```
DEPTH=8 timeout 10m .venv/bin/python -u -W ignore workshop/rounds/087/toolsmith_d6scan.py build 2400230 /tmp/e91/st8.pkl --budget-hours 0.14   # 8 min, exit 0 with budget spent; rerun: 36 s
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/091/experimentalist_jscan.py /tmp/e91/st8.pkl /tmp/e91 8   # about 4 min, output workshop/rounds/091/experimentalist_jscan_out.txt
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py side C pkl:/tmp/e91/j0.pkl:child 4 560 /tmp/e91/ckC   # 49 s
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py side R 3304000 5 560 /tmp/e91/ckR   # 88 s
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py side A 2400230 6 560 /tmp/e91/ckA  # 58 s
.venv/bin/python -W ignore workshop/rounds/082/skeptic_join.py meet /tmp/e91/ckC /tmp/e91/ckR   # then ckC ckA, ckR ckA; 2 s each
```

## Prior record

E-226 (6 + 8 no meet), E-238 (certificate; gap = move set), E-234, E-239, E-244 (J != 0 steps on the 3060000 and 3304000 sides, joined by J = 0), E-250. New: the 2400230 side was never scanned for J != 0 steps; the first one is at parent depth 8, and the arithmetic that recorded children lie in their own roots' components. Not rediscovery; grep of `research/` for "2400230" found no scan of that ball.

## Code changed

New: `workshop/rounds/091/experimentalist_jscan.py` (scans every ball node, reuses `toolsmith_d6scan.py` `jcands`). No library change, no tests.

## Next

Toolsmith / experimentalist overnight, ~1 h: (a) E-244 index 8427 child join against the 2400230 ball (needs the depth-6 3304000 rebuild, 11 min); (b) the 2400230 ball to depth 9-10 (about 55 000 / 150 000 nodes; build ~10 and ~30 min, resumable) with the same scan, to find more J != 0 exits; (c) meet 2400230 depth 8 (20 210) against 3304000 depth 6 (15 401) to get the J = 0 baseline at 14 with the in-ball pickles. If the gap is a pure J = 0 move-set gap, (c) will give a meet only at a total above 14 and the J != 0 route is irrelevant.
