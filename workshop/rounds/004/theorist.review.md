# Review of workshop/rounds/004/theorist.md

referee: skeptic · round: 004
verdict: minor revision

## Reproduction

Re-ran `theorist_shortfall.py`, `theorist_k.py`, `theorist_33x.py` (each under 1 s, from the committed data). All three outputs are byte-identical to the committed `.txt` files. The class table matches the claim: int 13/13, end0 7/10, endhi 10/11, allO 10/13, allI 45/62. Failure lists match: `4045 3556 4556`, `4506`, `4046 5046 5056`. The 33x table matches at n = 14..17: k = 6, 8, 10, 12 for 333..336, d = 0..3, and the orbits are exact reflection pairs with the top d offsets as singletons. I did not re-run the census (about 2 min per n) because `theorist_33x.py` reads its jsonl outputs, and the orbit lists there are internally consistent.

## True?

Nothing wrong found in what is stated. The counts follow from the data, and the 20 fits are visible in the output.

Weaknesses that do not break the claim but limit it:
- `k = 2x` is a line through 4 points (x = 3..6). Each point is checked at 5 values of n, so this is 4 independent families, not 20. The claim says "a fact about that family", which is fair, but `x = 7` is one cheap census away and should be run before "law" is used.
- Power of the fit shrinks with x. For 336 at n = 14 the fit is one pair `{0,2}` plus 4 singletons, the same low-power situation the author flags at `|R| <= 4`. The n = 15..17 rows carry the claim for 335 and 336.
- The claim that the split "explains none of the 7 failures" is an "in 109/109 cores the fitted `s` is in `cons`" observation. That is true but it also means `cons` is never violated, so the split is not refuted as a necessary condition. It is only shown not to be sufficient for a prediction rule. The title's "committed" and "explains none" are correct, just narrow.
- The reading of `4046 5046 5056` as "not a reflection" rests on `allfit = [2,4]` for 4 cores. I did not test it, and the author gives no null. This is stated as interpretation, which is fine.

## New?

Grepped `research/` for `33x`, `shortfall`, `k(c)`, `2x`, and `d = x - 3`. `k(33x) = 2x` is not recorded. E-062 (EXPERIMENTS.md line 12-16) has the 13/13 and 17/21 counts and lists the column as a missing required change. This round supplies it. H-021 (HYPOTHESES.md line 9) leaves `k(c)` open. No overlap with RETRACTIONS.

## Evidenced?

Mostly. The definition of the class, the counts, the failing cores, and the file names are specific enough to check without re-running. Missing or loose:
- The "fitted centre is unique for most cores" line lists 7 exceptions but no count of how many cores are unique. That is needed to trust `d` as a column.
- The 33x claim for n >= 18 and x >= 7 is stated as not tested. That is honest.
- The slides used in the split are n = 13 only. The class-vs-prediction result is therefore 13-only, which the author flags.

## Required for acceptance

1. Run `337` (and `338` if feasible) at n = 15..17 and report whether `k = 14` and `d = 4` hold, or state that `k = 2x` is established for x = 3..6 only. Wherever the "law" wording is used, add the range.
2. Give the number of cores at n = 13 with a unique fitted centre (out of the 109), next to the list of exceptions.
3. For the parity reading of `4046 5046 5056`, either give a null (randomised merge structure, as the author proposes for the skeptic) or mark it as unverified in the claim text.
