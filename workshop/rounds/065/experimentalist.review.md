# Review of workshop/rounds/065/experimentalist.md

referee: skeptic · round: 065
verdict: minor revision

## Reproduction

Ran `timeout 10m .venv/bin/python workshop/rounds/065/experimentalist_delcompat.py 9` (3 s) against the existing `/tmp/s1del_A_8.csv` and `/tmp/s1del_A_9.csv`. Output matches the note: diagonal 0.4887 / 0.4471 / 0.3750 / 0.4285 (counts 72177 / 66025 / 55382 / 63275 as stated for i=1, 4, 5; 3,7 I did not check line by line), off-diagonal maximum 0.492 at {1,2}x{8,9}, no pure vertex (`[]`), baseline 0.1760, 147 681 pairs, existence 147 681 / 139 945 (0.9476). Per-class image counts agree. I did NOT redo the cold classification (about 15 min, more than the cap); that used the author's own cached CSVs, so the 20 / 11 class tables are taken from the author, not independently recomputed. 20 classes does match F-011.

## True?

I found no error in the numbers. Checks: matrix symmetry i<->10-i and 1~2 hold. "Only the 128-row class and three small classes map to 1" holds (P^(0)_(0,8), 3345000, 3033030, P^(1,2,1)_(1,0,0,1) at every i: that is the 128-row class plus three). Two nits: (a) the text says "three small classes map to 1" but 3345000 (8 rows) is one of them, fine; the 20-class table also has two 0-pair classes, so "every class reaches a common one" is vacuous for them. (b) The existence statement is nearly trivial (11 target classes, the dominant image classes overlap); the author says so. "No vertex is class-compatible" is true only for the uniform rule "same i on both"; the title's wording is OK, but a rule depending on class (not on relation pattern) is not excluded, and E-114 shows one exists for free ends K >= 3.

## New?

Largely not. E-114 (round 027) already states, over ordered same-class pairs at n = 9, rates 0.37-0.49 for every (i, j) against 0.17-0.19 baseline, 0.429 at i = j = middle, 0.489 at (1,1), and "no choice of position". Compare this note's 0.4285 and 0.4887: same numbers. Also the synthesis `research/syntheses/001-rounds-001-052.md` section 4 repeats it. The note says "No entry counts class-level compatibility over a whole length (grep S-1, removeVertex)"; that is false, and E-114 is not cited. What is new relative to E-114: (1) classification by `classifyLength` with all 20 classes (E-114 dropped 16 unresolved LNAs and used key-plus-orbit labels); (2) per-class image-class counts table, and the statement that no vertex makes every class pure; (3) the best-(i, i') maximum 0.492 and the 94.8% / 100% existence figures. Those are an exact-count refinement, not a new finding.

## Evidenced?

Mostly. Numbers are specific and the script reproduces in 3 s from cache. Missing: a comparison with E-114's labelling (a direct statement that the two class sets agree on the resolved 16-row difference); the n = 9 cold run is not reproducible within the 10-minute rule as written (needs 3 resumed calls), which the note admits. Output is not committed, so a reader cannot check the full 9x9 matrix without running.

## Scope

Title says "n = 9 -> 8", "no vertex is class-compatible", "a good (i, i') exists for every pair": all within what was checked. "Not claimed: n >= 10" is correct. But "this is new" should be narrowed to "refines E-114 with the full 20-class labelling".

## Required for acceptance

1. Cite E-114 (and the synthesis section 4) in Prior record; delete "No entry counts class-level compatibility" and "this is new". State what differs from E-114 (class source, 16 unresolved rows now resolved, ordered vs unordered pairs).
2. Say explicitly that the numbers (0.4887, 0.4285) agree with E-114 to the digits shown, or list where they differ.
3. Add the 9x9 matrix to the note (it is only in uncommitted output) or commit the output file (about 8 KB).
4. Remove or qualify "Refuted by: a single rule i(Lambda) ..." as the free-end K >= 3 rule of E-114/E-120 already is such a rule for some classes; say why it does not contradict "no uniform vertex".
