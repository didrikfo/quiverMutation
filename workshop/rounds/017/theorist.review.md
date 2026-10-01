# Review of workshop/rounds/017/theorist.md

referee: skeptic · round: 017
verdict: accept (minor wording)

## Reproduction
- `theorist_rrule.py 14` (1 s): 250 present / 105 missing (all a=2), three-run 119/119. Matches.
- `theorist_class.py 13 34 35 36 37 55 235 455 2455 3334 444 2244` (6 s): sizes and J identical to the table (2386 {0,7}; 447 {1,6}; 449 {2,5} even o, 224 shadow at odd o; 501 {3,4}).
- Extra, not in the submission:
  - n=15: 3334, 2455 -> size 763 J={1,8}; 4445 -> 5648 J={0,9}. Matches the stated 763 / 5648.
  - n=17 (20 s, 4 s): 3334 -> 1191 J={1,10}, all 10 placements; 444 -> 11340 J={0,11}. Class 1 holds at n=17, and 3334 is not in class 0.
  - Fold: n=13 38 -> {3,4}, 39 -> {2,5}; n=15 38, 39 -> {4,5}. Consistent with x-4 folding by j <-> n-6-j at x=8,9. Also 3338 -> {2,5} (= 39) and 3335 at n=13 gives 449 {2,5} at odd o with a 224 shadow at even o.

## True?
No counterexample found. Lemma R is verified only against the rewrite engine over a bounded range (n=12, 14; a>=3; b<=7 or 8; d<=9), not proved. The author says "hand-derived, checked", and that is accurate. The a=2 exclusion is disclosed and sound.
Item 4 and the "class 1 != class 0" statement are explicitly scoped to this move set. Fine.
One small point: "x-4 at every placement (odd n)" is tested only up to x=9. I did not test x>=10 or letters >=6 in prefix words. The author does not claim these.

## New?
Grepped research/ for 3334, 2455, "444 -> 34", "lemma R", "class label". Only E-083/E-086/E-065/E-080 (memberships by size and by row set, and the drift label), as the author cites. The one-step reduction to 35, the shared lemma R and the J table are not recorded. Not in RETRACTIONS.

## Evidenced?
Yes. Ranges (n, a/b/d) and counts are stated, and the scripts and output files are named. Missing: no explicit statement that the J table was checked at n=17 (I did it, and it holds), and the "no invariant, closed under this move set only" caveat is correctly kept.

## Required for acceptance
1. Optionally add n=17 (3334 -> 1191 {1,10}; 444 -> 11340 {0,11}) to the evidence line.
2. State that R is verified by exhaustion over a bounded range, not proved, in the Claim heading ("proof" in kind overstates it).
