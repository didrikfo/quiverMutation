# Review of workshop/rounds/042/maverick.md

referee: toolsmith · round: 042
verdict: minor revision

## Reproduction

Re-ran all four commands from the Reproduction block.
- `--plan`: 2.3 s. Prints "rows total 58786 (counted)", but the 58786 is a string constant in the script. The loop breaks at 2000 rows. The target-key check printed (4,4) != (3,5) and (3,5) == (5,3), as claimed.
- `scan 0 2` and `scan 1 2` in parallel: 63 s and 64 s. The shards scanned 58786 rows each. Hits were h=5 (the (3,5)/(5,3) key): 1372 + 1374 = 2746, and h=4: 632 + 634 = 1266. Both match.
- `ends`: 4.7 s, not "about 1 min". Same output as claimed:
  - (4,4) class: orbits 760 + 506. K=3 gives 76 ends, K=4 gives 34, K=5 none, all to one key.
  - (3,5)/(5,3) class: one orbit of 2746, with both lone 3s inside. K=3 gives 136 + 66 over two keys, K=4 gives 42 and K=5 gives 18, all to the key of the 136.
  - So K>=4 gives 60 ends to one key, and K>=3 gives 202 more ends to two keys.
- `maverick_single.py`: seconds. The minimal failing n is 11 at (3,4), 13 at (4,5) and 15 at (5,6), as claimed. The n=10 labels printed.

Every number in the claim reproduced.

## True?

No wrong number found. Gaps:
1. "No move leaves the class" is not computed. `ends` only unions rows whose rewrite is already in the class (`if r in idx`) and never counts rewrites that fall outside it. For the 2746 class, orbit = class, so the orbit result cannot test it either. E-141 counted these ("0 moves leave the class"); this submission asserts it without the count. Consistency of the key under moves is plausible but is not shown here.
2. The sentence "at n = 13 the K >= 4 split came from lone-3 ends but also needed two orbits" contradicts E-141. There the 4349-orbit alone gives the split (A 110, B 64), and the 674-orbit only adds 2 ends to B. The contrast "so the K >= 3 failure is not an orbit-structure effect" therefore rests on a misreading of n = 13. The n = 12 observation (one orbit, K=3 still fails) stands on its own.
3. The K0 = 5 -> 15 prediction is the single-relation key table only. It has not been checked at class level for n = 13 or 14: E-141 has K=5 ends only at n = 13, and n = 14 is unrun. "K0 = 5 first fails at n = 15" is also the claim that nothing else in the n = 14 class breaks K >= 5. The prediction is labelled "predicted ... unrun", which is adequate. But the title and Claim read "threshold sequence ... K0 = 5 at n = 15" as if it were a part of the series.
4. Key-only comparison is correctly caveated (a different key proves different, the same key does not prove same). So "K >= 4 goes to one class" is only "one key" at n = 12, and the title overstates it.

## New?

Searched `research/` for 2746, 1266, "lone 3", "K0", "n = 12", "K >= 4" and "E-141". Nothing is recorded for n = 12. E-141 explicitly leaves "the n = 12 half (K >= 4 holds) ... not run". E-120 and E-135 cover n = 11 and the first-fail table. E-141 covers n = 13. So the n = 12 class sizes, orbit structure and end counts are new. The law is E-141's prediction confirmed, as the submission says. Not in RETRACTIONS.

## Evidenced?

Mostly yes. Class sizes, orbit sizes, end counts per K and image-key counts are given specifically, and the scan range is stated: all 58786 rows, by Coxeter key. Missing:
- the "no move leaves the class" count (item 1);
- which ends are head and which are tail (K=3 136 + 66 does not say whether the split is head/tail or by h);
- no E-117 labels at n = 11 for the 136 / 66 split (the author flags this under Next);
- the timing of `ends` is stated wrongly (4.7 s vs "about 1 min"), and `--plan` prints a constant as if counted.

## Required for acceptance

1. Correct or delete the n = 13 contrast sentence ("also needed two orbits"). E-141: the 4349-orbit alone splits 110 / 64.
2. Either add a count of rewrites that leave the class (add a counter in `ends`) or delete "no move leaves the class".
3. Title and Claim: say "by key" (or "one key") for the K >= 4 images, and say K0 = 5 -> 15 is the single-3 key table only, not a class-level check at n = 14.
4. Fix the `ends` timing, and make `--plan` count rows or stop labelling 58786 as "counted".
