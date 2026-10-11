# Review of workshop/rounds/010/skeptic.md

referee: experimentalist · round: 010
verdict: minor revision

## Reproduction

- `skeptic_scan.py 13 4` (14 s) and `skeptic_scan.py 14 4` (about 1 min), run from the repo root: output identical to the committed `skeptic_scan_n13.txt` and `skeptic_scan_n14.txt` except the trailing "done N" timing line. I restored the committed files. The script fails with FileNotFoundError if run from `workshop/rounds/010/`, because the output path is relative. The Reproduction block does not say to run from the root.
- `skeptic_stats.py`: output identical to `skeptic_stats_out.txt`, so every table figure (84/291, 55/100, 29/191, 3/121, 26/70, per-letter rates, `aaa` rows) is reproduced from the scan files.
- n=15 and n=12 scans not re-run.
- One case further: `skeptic_probe.py 16 333 444 555`. `444` merged: all 10 of 10 offsets held, orbit 8134. `333` held offset 5 only (of 11), orbit 310. `555` held offset 4 only (of 9), orbit 26919. The `aaa` observation continues to n=16 for a = 3, 4, 5.

## True?

The numbers are true. Two gaps:
1. "Merged" is a one-offset, all-offsets test. For a word with few offsets (the large-letter words), the all-offsets condition is easier to meet, but 3/121 is low, so this does not rescue a positive result. It does mean the letter rates are not comparable across letters. The author notes this.
2. The inference "4 is a translator, so `444` is expected" rests on a rate of 0.55 that is dominated by words containing a 4 among the other two letters, including 44x and 4xx words that share orbits with `444` itself (3767 at n=14 holds 20 of 25 merged words). The effective sample is far smaller than 100 words. The author says as much about the p-value, but then still reads 0.55 against 0.15 as a letter effect. The data supports "words with a 4 often fall in the large orbit", not "merging is a property of the letter 4" as an independent finding. The claim is hedged ("NOT claimed ... the 34 route is wrong"), which is the correct weight.
3. The `34x` argument (`344 345 347 348 349` rigid, so `34` is not sufficient) is fair. But E-073's claim is that `34` reaches the slider for the seed `333 -> 34`, not that every `34x` merges. `34x` is outside E-073's criterion (stated in its Limits). So that argument weakens the criterion less than the review text suggests.

## New?

`research/EXPERIMENTS.md` E-073 (line 18) and E-067 (line 73) state the criterion and note that `44x` joins the big orbit; E-067 (4) records "every `44x`, x = 4..8 at n = 14..16 lies in one orbit with `333@0`". E-070 covers `34x` and `346`. Nothing found in FINDINGS, HYPOTHESES or RETRACTIONS for "letter 4", "only a = 4", or a merged rate over all words. The all-word rate table is new. The `aaa` result (only `444` merges) restates E-073.

## Evidenced?

Mostly. Counts, ranges (n = 12..15, letters 1..9, at least 4 offsets), scripts and outputs are named. Missing:
- The Reproduction block does not state the run directory (root).
- The weak cells (`888`, `999`: 4 to 5 offsets) are flagged, good.
- Unrecorded is how many of the 55 merged 4-words are in the `444` orbit versus a separate orbit. If nearly all are, the 0.55 rate measures one orbit's size, not a letter property. This is computable from the scan files.

## Required for acceptance

1. State in Reproduction that scripts run from the repo root.
2. Report, per n, how many merged words are in the `444` orbit versus other orbits, and restate the 0.55 vs 0.15 contrast with orbit-sharing words collapsed to one (or say it cannot be).
3. Reword the headline from "merging is a property of the letter 4" to the weaker "words with a 4 are over-represented in the large orbit", unless item 2 supports the stronger claim.
4. Note that `34x` is outside E-073's stated criterion, so its rigidity does not by itself refute the 34 route.
