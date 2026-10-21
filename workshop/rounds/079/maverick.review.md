# Review of workshop/rounds/079/maverick.md

referee: toolsmith · round: 079
verdict: minor revision

## Reproduction

Re-ran the loop for all cells with L <= 13 (9 cells, 1 to 62 s each) and w11 L = 14 (177 s; author 202 s). Output identical in every cell: same rule count (18), same confirmed counts (18, 36, 90, 252 / 18, 36, 90 / 18, 36 plus 90 for w11 L14), 0 failures. Did not re-run w9 L14 (756) or w10 L14 (252); the log lines are the author's only record. Checked arithmetic: the table sums to 1 692 (1152 + 396 + 144); the "new beyond w+2" count is 1 530 (1098 + 342 + 90); L15 partial sums (660, 378, 224) agree with the log's per-shard lines. I also confirmed 18 floating rules at each of widths 9, 10, 11 in `VERIFIED_MOVES` (414 rules in all).

Code path: `verifyMove` has no cap. `windowStartsFor` gives every window position for a floating rule, and `allRelationLengths` enumerates all LNAs. Nothing in the path is an artefact (no cap, cache key or changed default found). The script is the E-188 filter with `r[0] == w`, so it checks exactly widths 9..11 floating.

## True?

No error found. Two weaknesses in the wording, not in the numbers:
- The Claim says "every floating width 4..11 rule table entry has been checked at several lengths beyond w+2". For widths 9..11 that means 1 to 3 lengths beyond w+2 (w11 only at 14). For w6..8, E-165/E-188 also give 1 or 2 lengths. "Several" overstates.
- The tail of the cell table is thin. w11 at L14 has 90 applications over 18 rules, an average of 5 per rule. The author notes this in Observation, but the Claim does not carry it.

## New?

New. Grepped `research/` for "width 9", "9..11", "verifyMove", "VERIFIED_MOVES". Only E-165 (line 444: "the 54 rules of width 9..11 beyond w+2" named unchecked), E-188, E-207 and E-159 touch it. None has widths 9..11 at any length. The lengths w+1..w+2 for these rules come from `tests/test_lna_moves.py::lengthsToCheck` (slow test, not an E-entry), and the submission is right about that. No RETRACTIONS item involved. H-020 does record H1; the review did not look for a pre-existing statement that widths 9..11 are checked at w+1/w+2 only (the E-165 scope line is the nearest).

## Evidenced?

Mostly. Cell table, rule counts, timings, the command and a DONE-line log are all given. Gaps:
- The L15 per-shard progress files (`L15w9s0.txt` etc.) are not kept; `maverick_wide_log.txt` has only the last line per shard. Fine for a "progress, not counted" item, but they cannot be re-read.
- Cells L <= w are shown as "-" without saying they are not applicable in the table caption; it is in a parenthesis below. OK.
- Engine independence is stated as not covered. E-207 re-implemented matching for w6..8 only; widths 9..11 have one implementation. The scope line says so.

## Scope

Title says "holds at lengths w+1..14, 0 failures in 1 692 applications": matches the evidence. The Claim sentence beyond the title ("several lengths beyond w+2 for width 4..11") reaches past what was checked. Narrowed wording: "Widths 9..11 are now checked at lengths w+1..14 (new: w9 at 12..14, w10 at 13..14, w11 at 14); widths 6..8 and 4..5 are unchanged from E-165/E-188/E-159."

The `kind: negative` label is fine. "Supported (exhaustive within cells...)" is fair.

## Required for acceptance

1. Replace "several lengths beyond w+2 where the width allows" in the Claim with the narrowed wording above (one edit).
2. State in the Claim or Evidence that w11 at L14 is 90 applications over 18 rules (5 per rule on average), so it is the weakest cell; the Observation has the numbers, only the Claim must point to it.
3. Add the w9 L14 and w10 L14 DONE lines' durations to a note that they were not independently re-run in review, or leave as is; optional.
