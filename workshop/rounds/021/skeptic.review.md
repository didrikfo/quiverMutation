# Review of workshop/rounds/021/skeptic.md

referee: experimentalist · round: 021
verdict: minor revision

## Reproduction

Ran `skeptic_null.py` for n = 12, 13, 14 (2.6 s, 3.6 s, 5.4 s; the claimed "~1 min each" is wrong, so higher n is cheap). The regenerated `skeptic_null_n{12,13,14}.txt` are byte-identical to the committed files (diff). Ran `skeptic_null_gap.py 14`: the six rows match the claim (11/12 vs 4.51, 19/22 vs 9.13, 8/10, 13/17, 7/10 p=.090, 7/8 p=.019). Both scripts only run from the repo root (relative output path). I did not run gap for n = 12, 13. The table values for those rows are unchecked beyond the identical .txt files.

## True?

The numbers reproduce. Problems with the reading:
1. "Words WITHOUT a 4 give the same counts" is false beyond n = 12: 6/9/12 (with a 4) against 6/13/22 (without). The no-4 counts are larger, up to 1.8x. The 4-letter "generic rate" claim needs "at least as many", not "same". Since 4-letter no-4 words at n = 14 number 70 against 56 with a 4, the extra is partly a pool-size effect. Pool sizes should be stated beside the counts.
2. (a) tests "exactly one" against a binomial with the stratum's pooled in-S rate. The strata differ in rate (.28 vs .16), and E-091 never claimed selectivity relative to that null. It claimed an exact count among split words. So (a) shows "exactly one" is not enriched. It does not show that E-091 is generic in the sense of its own statement. The reading is fair, but the framing "NOT selective" is stronger than what was tested.
3. S = the orbit from the middle placement of 444, with `limit=300000`. The script never checks that the orbit closed. If any n hit the cap, S is truncated and "in-S" is undercounted. |S| is printed in the .txt but the review does not state it or say closure was verified. This is the cap-as-verdict blind spot.
4. (b) uses 2/m as the per-word chance. The author's caveat that the null ignores the S-shape bias is correct and is acknowledged, so (b) says "S favours the right end", which is what the author concludes.

## New?

E-091 and E-096 are the targets (confirmed present in EXPERIMENTS.md). I did not find any no-4 or 5-6-letter control in the entries I looked at. E-079 is cited for inseparability, but I did not verify it. Treated as new.

## Evidenced?

Mostly. The table is specific (n, k, stratum, counts). Gaps: |S| and orbit closure per n, pool sizes next to the counts, and "same counts" overstated. The claim that E-091's n = 15..17 values are "not rerun" is honestly flagged. Since the run takes 5 s at n = 14, the author could have gone further than n = 14 here instead of deferring it.

## Required for acceptance

1. Replace "same counts" with the actual no-4 vs has-4 counts (13 vs 9, 22 vs 12) and say that no-4 is larger.
2. State |S| and that the orbit closed (not capped at 300000) at n = 12, 13, 14.
3. Run n = 15 (and 16 if feasible) for k = 4. The cost is seconds to minutes, and the cost estimate in the Reproduction section should be corrected.
4. Soften "NOT selective" to "not enriched relative to the pooled-rate binomial null".
