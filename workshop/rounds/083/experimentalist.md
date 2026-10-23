# S-1 sizing at class level: derivedOrbits(14) takes 319 s and 467 MB (fits), derivedOrbits(15) is predicted at about 20 min (would not fit); the key-level K-threshold law holds at n = 11, 12 and at n = 13 with K0 = 5 for (4,5) (E-141); orbit level undecided (n <= 13, lone-3 key classes)

author: experimentalist · round: 083 · kind: result
thread: S-1 · bears on: E-146, E-173, E-224
scope: lone-3 key classes only; key-level threshold per class = smallest K with one image key (n = 13 (4,5): 5, not 4); orbit level undecided, orbit != class (F-046); every n = 15 figure is a prediction, none measured; n = 11, 12, 13 laws; "class level" = `freeMoves.derivedOrbits(n-1)` orbit of the deleted image (a lower bound for the derived class: distinct orbits may be one class); ext = `derivedOrbits(..., edges=True, doubles=True)`, coarser, still only a lower bound; free ends are head/tail gaps K as in `maverick_n12.py`; n = 15 not run.

## Response to referee

Verdict minor revision. No re-run was needed: every correction below is a reading of the existing outputs `experimentalist_classlaw_{11,12,13}.txt` (the referee reproduced n = 11 and 12; the n = 13 run, 527 s, was not repeated, its file is unchanged).

1. n = 13 (4,5) threshold: done. The table has 2 keys at K>=4 and 1 at K>=5, so K0 = 5 (E-141: K>=4 transport fails). Claim, title, Prior record corrected; E-146 predicted 4 and E-141 refuted it. Also corrected from the same files: at n = 11, (3,4) has 2 keys at K>=3 and 1 at K>=4, so by the "smallest K with one key" reading that class gives 4, while (2,5) gives 3; the E-120/E-135 value K0 = 3 is a statement about the whole lone-3 key set, which I did not recompute. n = 12 (3,5): 4, as before. The old "threshold K" column was just the ends column; removed.
2. Done. "Does not reproduce" and "Refuted ... by" replaced by "orbit labels cannot confirm it; undecided" and "Decided by". F-046 (derivedOrbits orbits merge into classes at n = 11 without being moves) cited: orbit is not class. I add the referee's point that the plain orbit count falling with K (97, 17, 3, 3 at n = 12) is largely a coverage effect (fewer ends), not evidence of an orbit-level threshold; and that for n = 13 (3,6) the orbit count changes 15 -> 3 across K = 4 -> 5 while keys stay 1.
3. Done. From the n = 11 output at K>=3: (0,7) 1 orbit, (1,6) 1, (2,5) 2, (3,4) 3. Two of four classes, not three; the sentence is fixed. (My "three" counted K>=2 where (1,6) has 18 orbits.)
4. Done. All n = 15 numbers (1229 s, 1.7 GB, row count, 40 and 90 CPU-min, Stage A 6 min) now say "predicted"; Stage A pickle time is marked unmeasured and the plan total is stated without it. Sample check: n = 12/13/14 within 6%, the n = 14 miss (337 predicted vs 319 measured) is high-side, so n = 15 is an extrapolation one step beyond.
5. Not run this sitting (n = 11 `procedure` join test between the three n = 10 image orbits: needs a separate script and an unmeasured time). Stays in Next.
Also: E-224 sized the plan from a key-only n = 14 table; the Stage A here replaces that sizing, and says so below.

## Claim

Sizes (measured). `derivedOrbits(14)`: 742 900 rows, 188 093 orbits (91 144 singletons, largest 18 720), 319 s, peak RSS 467 MB, one process, so it fits one 10-minute command. `derivedOrbits(15)` (2 674 440 rows) was **not run**: the loop body (rewritesOf + stripLengthTwo) timed on 3000 evenly spaced rows per n predicts 21 / 81 / 337 s for n = 12 / 13 / 14 against measured 21.5 / 82 / 319 s, and (prediction, not measured) 1229 s (about 20 min) and about 1.7 GB (x3.6 factor) for n = 15. So the class-level n = 15 K0 = 5 scan needs a Stage A of one job that does not fit 10 min if it asks for `derivedOrbits(15)`, but it does not need that: only the orbit of the n = 14 images and the key class of the lone 3 are needed.

Law. At key level the E-146 K-threshold law holds at n = 11 ((3,4): K0 = 4; (2,5): 3), n = 12 ((3,5): K0 = 4) and n = 13 ((4,5): K0 = 5, (3,6): 4; E-141 refutes E-146's K0 = 4 for (4,5)): images of free ends with K >= K0 are one n-1 key, with K >= K0 - 1 two. At orbit level the images of the same ends fall in more than one `derivedOrbits(n-1)` orbit (table); orbit labels cannot confirm the law; undecided. Orbit is not class (F-046), so this is not a failure of the law for classes. Decided by: showing two image orbits are a single exact class (law then holds), or a K >= K0 image set splitting into two non-equivalent classes (law then fails at class level too).

## Evidence

Rows: the lone-3 key class = all rows with the Coxeter key of `0..3..0` (E-146 method); images by `removeVertex` at head and tail ends with K >= 2. Counts are distinct image labels over all such ends (ends / orbits / keys / ext orbits) for the class with the smallest (h,K) of the key class, i.e. the one carrying the threshold.

| n | key class (h,K) | ends | image keys | image orbits | ext orbits |
|---|---|---|---|---|---|
| 11 | (3,4) | K>=3 | 2 | 3 | 3 |
| 11 | (3,4) | K>=4 | 1 | 1 | 1 |
| 12 | (3,5) | K>=3 | 2 | 17 | not run |
| 12 | (3,5) | K>=4 | 1 | 3 | not run |
| 12 | (3,5) | K>=5 | 1 | 3 | not run |
| 13 | (4,5) | K>=4 | 2 | 7 | 3 |
| 13 | (4,5) | K>=5 | 1 | 2 | 2 |
| 13 | (3,6) | K>=4 | 1 | 15 | 1 |
| 13 | (3,6) | K>=5 | 1 | 3 | 1 |

Key level rows reproduce E-146 exactly (n = 12 class sizes 3837, 7057, 3378, 2746, 1266; 60 ends at K >= 4 for (3,5), 262 at K >= 3). Full tables: `experimentalist_classlaw_{11,12,13}.txt`. At n = 13 other key classes of the lone 3: (0,9) one orbit at every K; (1,8) one orbit from K >= 5; (2,7) 3 orbits, ext 2, from K >= 5. At n = 11, (2,5) has 2 orbits (ext 2) at all K >= 3 with one key. So at n = 11 (the cheapest case, 10 -> orbits plain 746, ext 113) two of four classes at K>=3 have more than one image orbit ((2,5): 2, (3,4): 3; (0,7) and (1,6): 1). Orbit counts falling with K (n = 12: 97, 17, 3, 3) are largely a coverage effect (fewer ends); at n = 13 (3,6) orbits go 15 -> 3 across K = 4 -> 5 while keys stay 1.

Shard plan (proposal only, nothing approved, nothing run; every n = 15 figure is a prediction): E-224's key-level plan reused. Stage A: one job, `derivedOrbits(14)` with orbit labels pickled (predicted about 6 min, 470 MB; pickle time unmeasured and included in no total). Stage B: 24 shards of the n = 15 enumeration, key per row only (E-224: 2.06 ms/row, about 3.8 min each), collecting rows whose key is the target lone-3 key class (5,6). Stage C: one job; for the collected rows (a few ten thousands predicted by scaling from 5023 at n = 13; unmeasured), form their move orbits (`maverick_n12.py ends` logic), delete K >= 5 and K >= 4 ends, look up the Stage A label. Stage C cost is small beside B. The per-row deletion keys of E-224's 62-min-per-shard variant are not needed at key level for this class, since only target rows are deleted. `derivedOrbits(15)` itself is a separate 20 min job and is not part of the plan.

Timing caution: the n = 13 key scan was 450 s alone (about 2.2 ms/row, above the 1.4 ms/row seen at n = 11); the key scan at n = 15 is the dominant cost and E-224's 2.06 ms/row is consistent with it.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/083/experimentalist_size.py 14         # 319 s (13: 82 s)
timeout 10m .venv/bin/python workshop/rounds/083/experimentalist_sample15.py        # about 90 s, predicts n = 15
timeout 10m .venv/bin/python workshop/rounds/083/experimentalist_classlaw.py 11 --ext   # 29 s
timeout 10m .venv/bin/python workshop/rounds/083/experimentalist_classlaw.py 12         # 120 s
timeout 10m .venv/bin/python workshop/rounds/083/experimentalist_classlaw.py 13 --ext   # 527 s (tight: 10 min cap)
```

## Prior record

E-146 (key-level law n = 12, prediction K0 = 4 at 13, 5 at 15), E-173 (n = 15 witness at key level, 33 min on 4 cores), E-224 (24-shard key plan), E-141 (n = 13 key level: K0 = 5 for (4,5), K>=4 fails; orbit-level split inside the 4349 orbit), F-046 (orbit is not class). Key-level n = 11, 12, 13 numbers are rediscoveries of E-146/E-141/E-120/E-135. E-224 sized the plan from a key-only n = 14 table; Stage A here (measured derivedOrbits(14)) replaces it. New: the measured cost of `derivedOrbits(13), (14)`, the n = 15 prediction, and the orbit-level counts above. Nothing in RETRACTIONS concerns these.

## Code changed

None in `quivermutation/`. New scripts in `workshop/rounds/083/`: `experimentalist_size.py`, `experimentalist_sample15.py`, `experimentalist_classlaw.py` (the `--save` flag writes a pickle to `/tmp/claude-0`, not the repo). No tests run (no library file touched).

## Next

Experimentalist/theorist: decide whether the 2 to 3 image orbits at the K0 level are one class (n = 11: three n = 10 orbits; exact `classifyLength(10)` is not available, E-206 used n = 9 labels). The n = 11 image orbits are the cheapest test: try `procedure` joins or the Smith-normal-form / Coxeter-polynomial hash between them. Toolsmith: pickle the Stage A table. Chair: n = 15 class-level scan is Stage A + 24 key shards of about 4 min, predicted about 40 CPU-min of Stage A/C (pickle time unmeasured, not included) and about 90 CPU-min of B, if approved. The n = 11 join test of the three image orbits is not run.
