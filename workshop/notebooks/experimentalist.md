# Experimentalist notebook (rewritten round 099)

## What I now believe
- 099 (T7): pair (8:4)(10:4) (starts differ by 2, overlap 2), k = 3 margin 2: 8 of 8 two-bystander cells stay at 2 (9-20 LNAs, 47-168 s each), BUT 0 controls lower: single bystanders (7,3), (12,3) share 2 arrows and stay at 2 at k = 3 and k = 4 (397-424 s). So "bystander sharing >= 2 arrows lowers" is false for this pair; the run of three needs consecutive starts. The cells are consistent with H-010 but uninformative until a lowering control exists.
- 097 (S-1 item 3): full n = 15 key scan (2 674 440 rows, ~46 min, four workers). Lone-3 class (5,6)/(6,5): 18 328 rows; free-end images at n = 14: 2 keys for K >= 3, 4, 5, 1 key for K >= 6. Class (4,7)/(7,4): 21 045 rows, 12 image keys: messy. Key level only.
- Row convention: n-vertex LNA is an n - 2 entry row. 094: pair (8:3)(9:3) k = 5 margin 1: 5 of 5 stay 2; margin 2 about 1 h (not approved). 086: k = 4 margin 2 cells 175-530 s.
- 091 (T10): 2400230 ball sizes 1, 9, 40, 130, 369, 1001, 2691, 7261, 20210; no bridge to 3304000. 089: |orbit(B_n)| = n 2^(n-5).

## What I tried
- 099: `rounds/099/experimentalist_pairs.py` (pair parameter), `_single.py`, `_single4.py`. Plan first with `--plan`.
- 097: `scan15`, `_batch.sh`, `_split.py`. Earlier: `k5`, `k5par`, `jscan`, `orb13`, `k4`, `_classlaw`, `_jkey`, `_join*`.
- Gotchas: shared machine (cores ~0.5 each); a 10 min timeout kills silently, a missing row is not a verdict; embedPattern returns None for overlapping patterns that do not embed (e.g. (9,3) next to (8:4)(10:4)); vacuous cells reach 0 LNAs; wait with an until-loop, not sleep.

## What I would do next
1. Find a lowering control for start-gap-2 pairs (stepped chain starts 8, 10, 12; other lengths) before any more two-bystander cells for them.
2. T7 m >= 3 at k = 4 or k = 5 margin 2: overnight only.
3. Move-orbit count of the 18 328 rows (097 next 1); S-1 item 1 at n = 10 -> 9.
