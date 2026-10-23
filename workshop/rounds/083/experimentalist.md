# S-1 sizing at class level: derivedOrbits(14) takes 319 s and 467 MB (fits), derivedOrbits(15) is about 20 min (does not); the K-threshold law holds at key level for n = 11, 12, 13 but only partly at orbit level (n <= 13, lone-3 key classes)

author: experimentalist · round: 083 · kind: result
thread: S-1 · bears on: E-146, E-173, E-224
scope: lone-3 key classes only; n = 11, 12, 13 laws; "class level" = `freeMoves.derivedOrbits(n-1)` orbit of the deleted image (a lower bound for the derived class: distinct orbits may be one class); ext = `derivedOrbits(..., edges=True, doubles=True)`, coarser, still only a lower bound; free ends are head/tail gaps K as in `maverick_n12.py`; n = 15 not run.

## Claim

Sizes (measured). `derivedOrbits(14)`: 742 900 rows, 188 093 orbits (91 144 singletons, largest 18 720), 319 s, peak RSS 467 MB, one process, so it fits one 10-minute command. `derivedOrbits(15)` (2 674 440 rows) was **not run**: the loop body (rewritesOf + stripLengthTwo) timed on 3000 evenly spaced rows per n predicts 21 / 81 / 337 s for n = 12 / 13 / 14 against measured 21.5 / 82 / 319 s, and 1229 s (about 20 min, about 1.7 GB by x3.6) for n = 15. So the class-level n = 15 K0 = 5 scan needs a Stage A of one job that does not fit 10 min if it asks for `derivedOrbits(15)`, but it does not need that: only the orbit of the n = 14 images and the key class of the lone 3 are needed.

Law. At key level the E-146 K-threshold law reproduces at n = 11 (K0 = 3), 12 (K0 = 4 for the (3,5) class) and 13 (K0 = 4 for the middle class (4,5)): the images of free ends with K >= K0 are one n-1 key, those with K >= K0 - 1 are two. At orbit level it does **not** reproduce as stated: the images of the same ends fall in more than one `derivedOrbits(n-1)` orbit (table), so the law cannot be read at this level; whether those orbits are one derived class is not decided. This does not claim the law fails for classes, only that orbit labels cannot confirm it. Refuted at orbit level by: showing two image orbits are a single exact class (it would then hold), or by a K >= K0 image set that splits into two non-equivalent classes (it would then fail at class level too).

## Evidence

Rows: the lone-3 key class = all rows with the Coxeter key of `0..3..0` (E-146 method); images by `removeVertex` at head and tail ends with K >= 2. Counts are distinct image labels over all such ends (ends / orbits / keys / ext orbits) for the class with the smallest (h,K) of the key class, i.e. the one carrying the threshold.

| n | key class (h,K) | ends | threshold K | image keys | image orbits | ext orbits |
|---|---|---|---|---|---|---|
| 11 | (3,4) | K>=3 | 3 | 2 | 3 | 3 |
| 11 | (3,4) | K>=4 | 4 | 1 | 1 | 1 |
| 12 | (3,5) | K>=3 | 3 | 2 | 17 | not run |
| 12 | (3,5) | K>=4 | 4 | 1 | 3 | not run |
| 12 | (3,5) | K>=5 | 5 | 1 | 3 | not run |
| 13 | (4,5) | K>=4 | 4 | 2 | 7 | 3 |
| 13 | (4,5) | K>=5 | 5 | 1 | 2 | 2 |
| 13 | (3,6) | K>=4 | 4 | 1 | 15 | 1 |
| 13 | (3,6) | K>=5 | 5 | 1 | 3 | 1 |

Key level rows reproduce E-146 exactly (n = 12 class sizes 3837, 7057, 3378, 2746, 1266; 60 ends at K >= 4 for (3,5), 262 at K >= 3). Full tables: `experimentalist_classlaw_{11,12,13}.txt`. At n = 13 other key classes of the lone 3: (0,9) one orbit at every K; (1,8) one orbit from K >= 5; (2,7) 3 orbits, ext 2, from K >= 5. At n = 11, (2,5) has 2 orbits (ext 2) at all K >= 3 with one key. So at n = 11 (the cheapest case, 10 -> orbits plain 746, ext 113) already three of four classes need more than orbit identity.

Shard plan (proposal only, nothing approved, nothing run): E-224's key-level plan reused. Stage A: one job, `derivedOrbits(14)` with orbit labels pickled (about 6 min, 470 MB; I did not time the pickle). Stage B: 24 shards of the n = 15 enumeration, key per row only (E-224: 2.06 ms/row, about 3.8 min each), collecting rows whose key is the target lone-3 key class (5,6). Stage C: one job; for the collected rows (a few ten thousands by scaling from 5023 at n = 13: unmeasured), form their move orbits (`maverick_n12.py ends` logic), delete K >= 5 and K >= 4 ends, look up the Stage A label. Stage C cost is small beside B. The per-row deletion keys of E-224's 62-min-per-shard variant are not needed at key level for this class, since only target rows are deleted. `derivedOrbits(15)` itself is a separate 20 min job and is not part of the plan.

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

E-146 (key-level law n = 12, prediction K0 = 4 at 13, 5 at 15), E-173 (n = 15 witness at key level, 33 min on 4 cores), E-224 (24-shard key plan), E-141 (orbit-level n = 13 split inside the 4349 orbit). Key-level n = 11, 12 and the n = 13 K0 = 4 numbers are rediscoveries of E-146/E-141. New: the measured cost of `derivedOrbits(13), (14)`, the n = 15 prediction, and the orbit-level counts above. Nothing in RETRACTIONS concerns these.

## Code changed

None in `quivermutation/`. New scripts in `workshop/rounds/083/`: `experimentalist_size.py`, `experimentalist_sample15.py`, `experimentalist_classlaw.py` (the `--save` flag writes a pickle to `/tmp/claude-0`, not the repo). No tests run (no library file touched).

## Next

Experimentalist/theorist: decide whether the 2 to 3 image orbits at the K0 level are one class (n = 11: three n = 10 orbits; exact `classifyLength(10)` is not available, E-206 used n = 9 labels). The n = 11 image orbits are the cheapest test: try `procedure` joins or the Smith-normal-form / Coxeter-polynomial hash between them. Toolsmith: pickle the Stage A table. Chair: n = 15 class-level scan is Stage A + 24 key shards of about 4 min, about 40 CPU-min of Stage A/C and about 90 CPU-min of B, if approved.
