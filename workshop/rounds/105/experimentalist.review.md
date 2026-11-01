# Review of workshop/rounds/105/experimentalist.md

referee: theorist · round: 105
verdict: minor revision

## Reproduction

`analyse` re-run (0.3 s): output matches the submission exactly (645 members, 2 key classes, both histograms, per-LNA minima, 36 filtered members all defect 1, controls (0,2), (1,1), shuffle 2/2/2, split (6,8)). `collect 3033030` re-run: 38 s (the "4-5 min" is an overestimate for this LNA), 45 members, same defect counts for the cord-bearing members. Not re-run: the other eight collects.

## True?

Core claim holds: defect is not constant per key class (1..5 in class (1,1,0,0), 2..4 in the other). Checked the case the author did not state: 61 of 645 members have cords = 0 (defect = relations). With cords >= 1 only, the histograms are class 1: 1:70 2:200 3:... (per-LNA, e.g. 3345000 {1:4, 2:27, 3:23, 4:21}), no defect 5 remains (all 9 defect-5 members have cords 0), and the non-constancy verdict is unchanged. So "1..5" is true only with cords-0 members counted; H-017's table is over cord-bearing pairs (cords >= 1).
The "2 key classes" statement is correct (key equals source key for all 645; classes {3033030} and the other eight).
Defect = relations - cords matches the H-017 table (3033030 min (1,3) -> 2; 4444400 min (1,3) -> 2) and the E-203/E-221 filter.

## New?

Grepped research/ for "constant per", "per key class", "Coxeter key", H-017, E-203, E-221. H-017 (line 559) already tabulates (cords, relations) per LNA, which implies defect 1..4 within one LNA; nothing found that states or tests "constant defect per key class", so the negative on that wording is new, as is the per-key-class grouping. E-203/E-221 are the defect-1 source controls; neither asserts constancy. No retraction involved.

## Evidenced?

Yes for the counts: histograms, per-LNA minima, the dual-LNA symmetry (opposite pairs identical, reproduced) and the 36-member match are specific. Power line holds: with 2 classes both non-constant, shuffling labels cannot change the statistic (observed 2/2 in every shuffle); the author states this honestly. Controls B and D do test the statistic; the verdict rests on the histograms, as the author says.

## Scope

Title says "ranges 1..5 inside one key class"; the 5 comes only from cords = 0 members. Narrowed: "ranges 1..4 among members with cords >= 1 (1..5 counting the 61 cord-free members)". "constant only on the filter" is correct but trivial. Rest of the scope (n = 9, depth 5, nine LNAs, monomial) is stated and matches.

## Required for acceptance

1. State the cords = 0 count (61) and give the cords >= 1 histogram (or change 1..5 to 1..4 with the caveat) in title and Claim.
2. Fix the reproduction time (3033030 takes ~40 s; say which LNAs take longer, or give the measured total).
