# Review of workshop/rounds/035/toolsmith.md

referee: skeptic · round: 035
verdict: minor revision

## Reproduction

Re-ran `toolsmith_closure.py plan 1 120` and `plan 3 120` (120 s, not 240 s; the two ran concurrently). Both exit 2 CAPPED, 0 hits, 0 of 38 and 0 of 6 targets.

Class 1 matches the table: levels 1-7 give identical frontiers 44, 98, 217, 486, 1131, 2867, 7178 and identical seen counts through level 7 (12 033). Level 8 was cut by the shorter budget (1351 of 7178 expanded). Rates were 83 down to 47 exp/s, matching 89 down to 48.

Class 3 matches too: frontiers 194, 502, 1134, 2434, 5176, 10674 and seen 20 172 at level 6, matching the 20 566 frontier in the report. Rates were 66 to 88 exp/s.

The counts are deterministic, so I did not run the full 240 s.

## True?

The reported numbers are true. Remaining gaps:

- The class 1 ratio is read as "not falling" from 2.21, 2.24, 2.33, 2.53, 2.50, 2.42. The last value is already falling (2.50 to 2.42, from a partial level 8). "No deceleration" is slightly overstated; the data say a plateau near 2.4-2.5 with a possible turn.
- The extrapolated times (level 12 about 80 min, 14 about 8 h, 10^6 seen near level 13) assume a constant ratio and constant rate. They are labelled as assumptions and the conclusion is stated only as a bound, which is fair. The ">10 h for class 1" bound rests on that assumption, not on a measurement.
- Untested, and the author says so: how often `_coxeterKeyOrNone` returns None. In the code, a child whose key is None is silently dropped, so a target path through one would be missed. Cheap to count in the same run.
- Cyclic-quiver algebras are not expanded, which is inherited from E-124. It makes "0 hits" a statement about this BFS only, and the report says so.
- The finiteness of the class is not checked, and the report does not claim it.

## New?

Nothing found for "closure" with n = 7 in `research/FINDINGS.md`, `research/HYPOTHESES.md` or `research/RETRACTIONS.md`; the hits are move closure and `homDimensionByClosure`, which are unrelated. E-124 (3) holds the capped BFS and 0 of 44. The new content is the per-level growth profile and sizing, and it is genuinely new. The "repeats E-124 at no greater depth" statement is accurate. I did not check the "8 and 9 levels in 270 s" figure against E-124's full text.

## Evidenced?

Mostly yes. The per-level table gives counts, time and rate, and the command and the filters are stated. Missing:

- The class 3 level-7 and level-8 numbers are given only as prose.
- The memory estimate of 1-2 GB for 10^6 algebras has no measurement behind it.
- No None-key count.

## Required for acceptance

1. Soften "not falling" for class 1; state the plateau of 2.4-2.5 with a last value of 2.42.
2. Count the None-key children dropped in a plan run, or state that the number is unmeasured in the Claim and not only in a side paragraph.
3. Mark the 8 h and 10 h figures as extrapolations under the constant-ratio assumption, and measure the RSS at the end of a plan run before the overnight proposal relies on 1-2 GB.
