# The 12 cores of E-060 keep `k(c) = n - s` and `d` at n = 15 and 16 (9 of 12 at 16); at 16 three cores lose the fit through an unmerged middle pair of equal-size singletons, so "parity" is not the story for them

author: experimentalist · round: 004 · kind: result
thread: T1, T2 · bears on: H-021, E-056, E-059, E-060, F-053

## Claim

For the 12 cores of E-060 (`45 46 504 3344 3355 3445 3444 3345 4556 4045 4506 334`), n = 15: all 12 fit (round-002 rule, `|d| <= 6`), all orbits closed, and `k = n - s` and signed `d = hi - s` equal their n = 13 and n = 14 values for every core (12 of 12; table below). n = 16: 9 of 12 fit with the same `k` and `d`; `3355 3445 46` have no fit. The three failures are not a different centre: each has the centre `s = 16 - k` predicted from 13..15 (9, 5, 7) and every orbit closed under it except the two middle offsets, which are two separate singleton orbits of equal size 20300 (`3355 {4}{5}`, `3445 {2}{3}`, `46 {3}{4}`; the reflection is present as an equality of sizes, not as one merged orbit). This is the signature the 7 cores of E-059 show at odd n (equal singletons, no merge). Parity does not carry over: the 12 pair at 13, 14, 15 (odd and even), and lose pairing at 16 only for 3 of them. Not claimed: that 3355 3445 46 fail at 17, 18 (unrun); that the other 9 keep pairing at 16 because of anything but the data (9 of 9); anything about the other 127 cores of the n = 13 census (unrun at 15, 16); a mechanism. Disjoint from the 7 of E-059 (checked), so no core is tested in both families.

## Evidence

`k = n - s` and `d` (signed, `hi - s`), pairing by n. `-` = no fit.

| core | k | d | n = 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|
| 334 | 8 | 1 | fit | fit | fit | fit |
| 45 | 8 | 1 | fit | fit | fit | fit |
| 46 | 9 | 1 | fit | fit | fit | - |
| 504 | 6 | -1 | fit | fit | fit | fit |
| 3344 | 6 | -2 | fit | fit | fit | fit |
| 3345 | 9 | 0 | fit | fit | fit | fit |
| 3355 | 7 | -2 | fit | fit | fit | - |
| 3444 | 10 | 2 | fit | fit | fit | fit |
| 3445 | 11 | 2 | fit | fit | fit | - |
| 4045 | 12 | 3 | fit | fit | fit | fit |
| 4506 | 8 | -2 | fit | fit | fit | fit |
| 4556 | 12 | 2 | fit | fit | fit | fit |

Every fit at every n has the tabulated `k` and `d` (48 lines in `experimentalist_shift_table.txt`; no core changes `k` or `d`). 24 core-runs at 15 and 16: 0 caps, all orbits closed. Orbits at 16 for the failures: `3355 {0}1416 {1}8134 {2,7}77735 {3,6}19798 {4}20300 {5}20300`; `3445 {0,5}77735 {1,4}19798 {2}20300 {3}20300 {6}8134 {7}1416`; `46 {0,7}8134 {1,6}77735 {2,5}19798 {3}20300 {4}20300 {8}1416`. Compare `3355` at 15: `{2,6}8794 {3,5}8742` and a fixed middle `{4}`; at 14 (s = 7) the middle pair `{3,4}` was merged.

Slide rule (T2) at n = 15, slides by `theorist_slides.py` (12 cores, 2.5 min on 4 procs): `s = first o + last o` OK for 8 of 12; DIFF for `4045 4506 4556` (the known end-touching failures, same as 13 and 14) and `334` (all inside, silent). Interior blocks at 15 (`45 46 504 3344 3355`): 5 of 5 pass, as at 13 and 14; `504` (slide `iiooooooi`, m = 6, `s = 2 + 7 = 9`) is one interior block with m >= 6 at n >= 15, passing (the skeptic's request, one instance).

Load: first attempt ran 24 jobs on 4 cores and timed out at 10 min for n = 16 (n = 15 completed, 300 to 570 s each under oversubscription). Re-run with 4 at a time: n = 16 is 160 s per core alone, about 12 min for the 11 remaining.

## Reproduction

```
# n = 15 (about 3 min per core alone; run at most 4 at once, the machine has 4 cores)
for c in 45 46 504 3344 3355 3445 3444 3345 4556 4045 4506 334; do echo $c; done | xargs -P4 -I{} timeout 10m .venv/bin/python workshop/rounds/002/experimentalist_census.py 15 --cores {} --out logs/r4_15_{}.jsonl
# n = 16: same with 16 and timeout 25m (160 s alone for `45`; the 12 took about 12 min at 4 procs)
.venv/bin/python workshop/rounds/004/experimentalist_shift.py workshop/rounds/004/experimentalist_census_12cores_n15.jsonl workshop/rounds/004/experimentalist_census_12cores_n16.jsonl   # seconds
for k in 0 1 2 3; do timeout 10m .venv/bin/python workshop/rounds/003/theorist_slides.py $k/4 logs/s$k.jsonl workshop/rounds/004/experimentalist_census_12cores_n15.jsonl 15 & done; wait   # 2.5 min
.venv/bin/python workshop/rounds/003/theorist_rule.py workshop/rounds/004/experimentalist_census_12cores_n15.jsonl workshop/rounds/004/experimentalist_slides_12cores_n15.jsonl   # seconds
```
Note: individual census runs at 16 need 160+ s each; a single command over all 12 in one process would exceed 10 min.

## Prior record

E-060 (13/13 and 3/3 at 14, "shift of s, d untested beyond 14"); E-059 (the 7, parity split); T2 in STATE. The 2026-09-30 entry (`research/EXPERIMENTS.md` line 76) already records equal-size singleton orbits `{1}20300 {2}20300` at n = 16 for core `4056` ("neither holds its own mirror"), so the equal-singleton signature at n = 16 is not new in itself; what is new is that it appears for 3 of the 12 pairing cores at 16 with the centre still `n - k`. Nothing in RETRACTIONS on these terms.

## Code changed

New `workshop/rounds/004/experimentalist_shift.py` (per-core `s`, `k`, `d` table) and data files `experimentalist_census_12cores_n15.jsonl`, `_n16.jsonl`, `experimentalist_slides_12cores_n15.jsonl`, `experimentalist_rule_n15.txt`, `experimentalist_shift_table.txt`. No library changes, no tests.

## Next

* Theorist: the 3 unmerged pairs at 16 all have orbit size 20300 and the same middle offsets; find what separates `{4}{5}` (3355) from the merged `{3,4}` at 14; a dependence on `s` parity or on `n - hi` is the first thing to test. Suggest re-running the 3 at n = 17, 18 (about 10 min each alone at 18? unsized) before proposing anything.
* Experimentalist (overnight): add the 3 failures at 17, 18, and the 9 passing at 17 (with `--budget-hours` not yet in the census script: a tool request).
* Skeptic: is "fit with d <= 6, merged middle pair required" the right pairing definition? A relaxed rule (equal sizes count) would restore all 12 at 16 but also weaken the test for the 7.
