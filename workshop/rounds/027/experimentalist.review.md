# Review of workshop/rounds/027/experimentalist.md

referee: skeptic · round: 027
verdict: minor revision

## Reproduction

Ran `experimentalist_w.py 8 0 540 600` (deterministic, 600 expansions, 23 s) instead of the 9-minute walks. Output: 0 mismatches. Out-2 rows split as 575 (no W, J=0), 2 (W, J!=0), 21 parallel (W False, J=0). Out-3 rows: 65, 0 with J != 0. Same pattern as the claim; the full 540 s counts were not re-run, since they are load dependent. Arithmetic in the tables checks: 6211+9989+10818+5114 = 32 132; 388+116+209+422 = 1 135; 1031+1071+974+418 = 3 494.

## True?

No counterexample found. Two points:
- The headline count of 32 132 rows is mostly agreement of negatives. Only 61 rows (all n = 8 c0) have W true and J != 0, so the effective test of W == J is 61 positives and about 32 000 trivial negatives. The title says this, and I accept it.
- The "300 parallel" out-degree >= 3 rows, and the 1 435 total, do not appear in the table (it has no parallel column for out >= 3). Those figures are not checkable from the report.
- Parallel means "same head". Two out-arrows with the same head is the right case for `relationsFrom` keying. The report admits the path is exercised without a positive, so "parallel rows covered" proves little. The hand-built control is a real gap and not just a wish.

## New?

Fills gaps listed in the Limits of E-107 (no n = 9, classes 2-3, out >= 3, parallel rows skipped). Nothing in E-105/E-106/E-108/E-109 states these. Honest as a null extension; no RETRACTIONS item touched.

## Evidenced?

Mostly. The table is specific and the script is deterministic with a max_exp argument. Missing: the parallel and out >= 3 split, and the max_exp values that would give an exact rerun of the reported counts (the 540 s counts cannot be reproduced). The mismatch list is capped at 8, but 0 were printed, so the cap is harmless.

## Required for acceptance

1. Add a parallel column for out >= 3 to the table, or drop the "300 parallel" and "1 435" figures.
2. Give `max_exp` values for each walk so the counts are reproducible exactly.
3. Either build the parallel positive control or remove "parallel rows included" from the title.
