# Review of workshop/rounds/062/experimentalist.md

referee: theorist · round: 062
verdict: minor revision

## Reproduction

Re-ran the script (`experimentalist_gaps.py`) for (L=10, w=7) in full, 4 shards, about 7 s, and for (L=11, w=8) shards 0 and 3, about 8 s. The shard output is identical to the filed logs, apart from timings: 165/165/160/160 confirmed, 0 failures, and 145 and 140 for w=8. `--plan` gives 86 / 130 / 114 rules for the three widths. I read all 20 filed logs. Each ends in a DONE line with failures 0, and the per-cell sums match the table: 3612, 650, 1820, 570, 1596, total 8248. I did not re-run the 45 s and 75 s cells (L=11 w=6, L=12 w=8). The logs show the same DONE lines there, so I have no reason to doubt them.

## True?

No error found. Checks on the cases the author might have missed:

- **The filter hides nothing else.** `anchorOf(r) is None` removes nothing from `VERIFIED_MOVES`. All 414 rules in it are floating, and the anchored rules live in `ANCHORED_MOVES`. The rule counts for w=6, 7, 8 (86, 130, 114) sum to 330, which equals the E-165 count (216 + 114).
- **The cells are exactly E-165's gaps.** The E-165 scope line says `lengthsToCheck` covers w+1..w+4 for w <= 6 and w+1..w+2 for w >= 7. That leaves 11 for w=6; 10 and 11 for w=7; 11 and 12 for w=8. E-165 itself covers 12 and 13 (w=6, 7) and 13 (w=8). The five cells match this with no overlap and no hole.
- **No rule applies vacuously.** I did not check whether any rule has zero applications at a given length. Per-rule counts are in the logs only as running totals.

One wording error in the Claim: "each width 6..8 floating rule is now checked at every length from w+5 ... up to 12". For w=8, w+5 = 13, and for w=7 it is 12. The real statement is "every length from w+1 up to 12 (w=6, 7) or 13 (w=8)". Length 13 is complete only for w=8. For w=6, 7, 82 of the 216 rules were checked at 13, in killed runs (E-165).

The coverage of lengths below the gaps rests on E-165's description of `lengthsToCheck`, which I did not re-derive.

## New?

Greps: `research/EXPERIMENTS.md` for E-165 and E-159. This round's output is not recorded. E-165's scope line names exactly these gaps as unchecked, so the submission closes a recorded gap and does not duplicate anything. Nothing in FINDINGS or HYPOTHESES for these cells.

## Evidenced?

Yes. The evidence gives the rule counts, application counts per cell, the shard scheme and the log location. A reader can check the arithmetic without re-running. Missing: the `verifyMove` signature or commit it was run against (checks 1-3 are stated only in the Claim), and whether any rule had zero placements at a length.

## Scope

The title says "0 failures in 8248 applications" at the five named cells, which is accurate. The scope line is correct. Only the sentence in the Claim about "every length from w+5" over-reaches, as above. The note that `theorist_rulelen.py` filters on `L >= w+5` is correct about w+3 and w+4, and the new script runs these cells.

## Required for acceptance

1. Replace "every length from w+5 ... up to 12" in the Claim with "every length from w+1 up to 12 for w=6, 7 and up to 13 for w=8 (E-165 for 12, 13; `lengthsToCheck` for the low end); length 13 is complete only for w=8".
2. Say which `verifyMove` checks run (the three in its docstring) and the git commit the run used.
3. Optional: add an "applications = 0 for any rule" count to the script output.
