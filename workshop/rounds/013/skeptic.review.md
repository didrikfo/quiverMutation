# Review of workshop/rounds/013/skeptic.md

referee: experimentalist · round: 013
verdict: minor revision

## Reproduction

Re-ran `skeptic_orbscan.py` for n = 12, 13, 14, 15 in parallel (output redirected to scratch, author's files untouched), about 4 minutes in all. `diff` against the committed `skeptic_orbscan_n12..15.txt`: identical for all four (35/56/84/120 lines). `skeptic_orbstats_out.txt` read against the scans: every per-n count, the pooled table, the big-orbit lists and "7/8, 9/12, 11/16, 13/19" match. I checked by hand n = 14: merged words are in orbit 3767 (11 words), 886 (`236 266 466`), 491 (`238 288 488`), 11820 (`458 468 568`) and single orbits `233 255 277 457`, as claimed. `344 345 347 348 349` rigid, `346`, `234` merged in the big orbit at n = 14, 15.

## True?

Yes as to the data. One discrepancy with the record, not with the scan: E-077's Limits say "at n = 14 the orbit 3767 holds 20 of the 25 merged words". The scan has 11 of 25 merged words (24 without `222`) in 3767. Either E-077 counted something else (offsets? held counts?) or it is wrong; the submission should say which, since it claims to fill that gap. Unchecked: `aaa` for a = 3..9 only in the 4-n pool; the "orbit" is the orbit at the middle offset only, so "orbits hosting a word" for rigid words is a convention (a rigid word spans several orbits). The `{4aa, 2aa}` orbit-mate pattern is listed, not tested against a rule. No n = 16 and no 4-letter words (stated).

## New?

Partly. E-077 (Limits) already says the 0.55 rate "measures largely one big orbit" and asked for the count of merged words in the `444` orbit; this supplies it, plus `234`/`346` in the `444` orbit and the 2-driven small orbit list. Grepped `444`, `34x`, `orbit 3767`, `E-077`, `letter 4` in EXPERIMENTS/FINDINGS/HYPOTHESES: nothing else. No RETRACTIONS entry. "Refutes the 55/100 vs 3/121 as an effect size" is stronger than E-077's own caveat; it is a sharpening, not a refutation.

## Evidenced?

Mostly. Specific (n range, word class, offset threshold, drop of size-1 cells, file names, counts). Weaknesses: (a) the pooled orbit table counts the same big orbit once per n and the author says so, so orbit rows are not independent either; no per-n orbit-level contrast is given for no-4 vs 4 beyond counts. (b) Claim 2 (cannot separate letter-4 from collapse-to-`34`) is a statement about this design, correctly hedged. (c) "kind: negative / refutes" overstated for item 4.

## Required for acceptance

1. Reconcile "20 of 25" (E-077) with 11 of 25 for orbit 3767 at n = 14, and state which is right.
2. Change "refutes" to "sharpens" in the title/header, or justify.
3. Say in item 4 that "orbits hosting" uses the middle offset's orbit only.
