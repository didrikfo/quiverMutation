# Review of workshop/rounds/003/experimentalist.md

referee: theorist · round: 003
verdict: minor revision

## Reproduction

* `experimentalist_table.py` on the saved 74-row JSONL: 0.3 s, output matches the submitted table row for row (n = 8..18).
* Fresh census, not from the saved file: cores 344 and 4405 at n = 14, 15, 16 (`experimentalist_census.py`, 3.4 to 4.2 s each). Orbits are identical to the saved rows. `experimentalist_fit.py` gives 344: d = 0 s = 7 (n = 14), no fit (n = 15), d = 0 s = 9 (n = 16); 4405: d = 2 (n = 14), no fit (n = 15), d = 2 (n = 16). The quoted orbit lists for 344 (`{2,5}114` at n = 14; `{2}64 {6}64` at n = 15) and the `{0,s}` sizes 47, 52, 57 reproduce exactly.
* Not re-run: the other 5 cores, and n = 8..13, 17, 18 (the saved file is consistent with the runs I did, but I did not regenerate it). The n = 13 rows are from the round 002 census, as the author says.

## True?

Nothing found wrong in the range stated (7 cores, n = 12..18). Two limits:

* Selection: the 7 cores were picked because they show strict mirror without pairing at n = 13. "Pairs at 12, 14, 16, 18 and not at 13, 15, 17" is 7 of 7 for those cores, which is a real regularity, but it says nothing about how the other 132 cores split by parity. The title's "the defect is a parity effect" is a claim about these 7 only; the claim paragraph states this, the title does not.
* "No fit" at odd n is a statement about the fit rule (closure under o -> s - o with s within 6 of lo + hi), not about the walks. The gloss that the odd-n singletons "would pair" (equal sizes 50/50, 64/64) is interpretation; equal orbit size is not a reflection. The evidence for it is that the same equalities appear at 15 (64/64) and, for 4405 at n = 15, `{0}64 {4}64`, `{1,3}51`, `{2}68`, but no mechanism links those.

## New?

grep of `research/` for parity, odd n, even n, odd length, even length: nothing relevant (one unrelated hit in FINDINGS.md line 229 about arrows). Related and not duplicates: E-056 (the 7 at n = 13, "pairing with a defect"), E-052 and H-021 ("eight cores pair at n = 13 and 14"; this submission shows that for these 7 the n = 13 non-pairing is the odd exception, and extends the pairing to 12, 16, 18). It also bears on F-053. Nothing in RETRACTIONS. The result is new.

## Evidenced?

Yes for what is stated: the table gives n, cores, pair, strict, no-fit; the orbit listings are specific; the "not claimed" list is honest (n > 18, n < 12, other 132 cores, mechanism). Gaps:

* The 74-row file and the loop leave the 8-core versus 7-core discrepancy unexplained: H-021 says eight cores pair at n = 13 and 14, E-056 says 8 cores have strict mirror true. The submission does not say which E-056 strict-mirror core is not among the 7 (or whether 8 = 7 + `3346`-like).
* d = 0 for 4 of 7 at n = 14 is the strong evidence; d = 1, 2 for the rest rests on the generous d <= 6 rule. The author states this and passes it to the skeptic; fine, but the table should give d per core per even n (a column, or a saved text file), since the "7 of 7 pair" figure is only as strong as the d values.
* Title should read "for these 7 cores".

## Required for acceptance

1. Retitle or qualify the title: the parity effect is shown for these 7 selected cores, not for cores in general.
2. Add d (and s) per core at each even n to the evidence (a small table or a saved text file next to the JSONL), so "pairs" is checkable at a glance.
3. Say how the 7 relate to the 8 cores of H-021/E-052 and the 8 strict-true cores of E-056 (which one is missing and why).
4. Mark the "would pair" reading of equal-size singletons as interpretation.
