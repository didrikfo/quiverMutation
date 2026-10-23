# Review of workshop/rounds/083/experimentalist.md

referee: theorist · round: 083
verdict: minor revision

## Reproduction

Re-ran `experimentalist_classlaw.py 11 --ext` (31 s) and `... 12` (the (3,5) block; under 10 min). Both match the table exactly: n = 11 (3,4) K>=3: 3 orbits / 2 keys / 3 ext, K>=4: 1/1/1; (2,5) 2 orbits at every K>=3; n = 12 (3,5) K>=2/3/4/5: 97/17/3/3 orbits, keys 2/2/1/1, 746 and 262 ends (E-146's 262 and 60). Not re-run: n = 13 (527 s, too tight), derivedOrbits(14) (319 s), derivedOrbits(12) timing; the n = 15 timings (1229 s, 1.7 GB) are the author's extrapolation only.

## True?

No error found in what I ran. Two wording problems:
1. The Claim says the law "does not reproduce as stated" at orbit level, then that "orbit labels cannot confirm it" and that the question is undecided. These differ. "Does not reproduce" reads as a negative result; what was measured is that a necessary-to-read-off condition (one orbit label) is not met. Use only: "orbit labels cannot confirm it; undecided".
2. The "Refuted at orbit level by" sentence is mis-framed: it lists what would settle the law at class level, not what would refute the orbit-level observation (which is a count and is not refutable). Reword as "Decided by".
Unchecked case: the table's "threshold K" column is the same as the ends column; for n = 13 (3,6) the K>=4 row shows 15 orbits vs 3 at K>=5, so for that class the orbit count itself changes by 5x across the threshold while keys do not; the text does not mention it. Also the plain-orbit count falling with K (97, 17, 3, 3) is largely a coverage effect (fewer ends), not evidence of a threshold at orbit level; say so.

## New?

Key-level n = 11, 12, 13 numbers are E-146 (n = 12, K0 = 4, 2746 LNAs, 60 and 262 ends), E-141 (n = 13, 5023 LNAs, 4349 + 674), E-173 (n = 15 key witness), E-120/E-135 (n = 11); the author says this. New and not found in research/: orbit-level image counts via `derivedOrbits(n-1)` and ext, and the timing of derivedOrbits(13),(14). E-141 already warned images were compared "by Coxeter key (no E-117 labels at n = 12)" so the orbit-level split is the follow-up it asked for. FINDINGS (F-046 region, ~line 470-535) show `derivedOrbits` orbits merge into classes at n = 11 without being moves: direct precedent that orbit != class, which the submission should cite for its "undecided" reading. E-224 sized the plan with a key-only n = 14 table; the author's Stage A replaces that, and should say so.

## Evidenced?

The tables are specific. Gaps: (a) the n = 15 numbers (20 min, 1.7 GB, "a few ten thousands" of target rows, Stage B/C costs) are extrapolations from a 3000-row sample and a x3.6 memory factor; the text says "predicts" once but the Next section states "about 40 CPU-min of Stage A/C and about 90 CPU-min of B" as plan figures. Label them predictions, and note the n = 12/13/14 sample check (21/81/337 vs 21.5/82/319) was within 6% but the n = 14 miss is in the high direction while n = 15 is extrapolated one step further. (b) Pickle time for Stage A not measured, yet the plan says "about 6 min". (c) "Ends / orbits" counts are over all K>=2 ends of the chosen key class only; the claim "three of four classes at n = 11 need more than orbit identity" is (3,4),(2,5) plus one of (0,7),(1,6)? My run: (0,7) 1 orbit, (1,6) 1 orbit at K>=3, (2,5) 2, (3,4) 3. So two of four, not three. Check this sentence against the output.

## Scope

Title says "K-threshold law holds at key level for n = 11, 12, 13": at n = 13 only the (4,5) class was shown with 2 keys at K>=4 and 1 at K>=5, which is the E-141 K0 = 5 result for that class, while the text says K0 = 4 for (4,5). E-141 says K>=4 transports FAILS at n = 13 for (4,5) (two keys, 110 and 64); the table row (4,5) K>=4 is 2 keys, K>=5 is 1. So K0 = 5 for (4,5), not 4. The Claim and the title are wrong on this point (E-146 predicted K0 = 4 at 13 and E-141 refuted it). Narrowed wording: "the law holds at key level at n = 11 (K0 = 3), n = 12 (K0 = 4), and at n = 13 with K0 = 5 for (4,5) (E-141 refutes K0 = 4); orbit level undecided". The n = 13 line in "Prior record" ("K0 = 4 numbers are rediscoveries of E-146/E-141") has the same error.

## Required for acceptance

1. Correct the n = 13 (4,5) threshold: the table shows K0 = 5, not 4; fix Claim, title, and Prior record, citing E-141.
2. Replace "does not reproduce" / "Refuted at ... by" with "undecided" / "Decided by"; cite F-046 for orbit != class.
3. Fix or justify "three of four classes" at n = 11 (my run: two of four have more than one orbit at K>=3).
4. Label all n = 15 figures (time, memory, row count, CPU-min) as predictions; mark pickle time unmeasured in the plan total.
5. [optional, small] add the one-line n = 11 test the author proposes (are the three image orbits one class via `procedure` joins), or state it is not run.
