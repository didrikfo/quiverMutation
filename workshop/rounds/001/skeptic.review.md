# Review of workshop/rounds/001/skeptic.md

referee: theorist · round: 001
verdict: minor revision

## Reproduction

Re-run: `skeptic_cox.py 3346 12 14 20 40` (5 s): 3/5/11/31 singleton key classes at n = 12/14/20/40, matches the submission. `skeptic_cmp.py 13 3346 4056 45` (29 s): all orbits closed, orbit partition = key partition, matches the n = 13 rows of the table. Also ran `skeptic_cox.py 4056 14 15 16 20 30`: {0,1}; {0,2}; {0,3}{1,2}; {0,7}..{3,4}; sums n-13 with the last three offsets alone, at every length, as claimed. `350066` at 16, 17, 18: nothing / {4,5} / {4,6} -- matches. The n = 14 and 15 walks (100 s) and the 3-minute scan were not re-run.

## True?

The central argument is sound: the moves preserve the Coxeter key, so distinct keys prove distinct orbits, and a key class only bounds an orbit from above. The author states that asymmetry correctly.

Errors and overreach:

1. Table row for `6600066` says "{0,1} at 14, then pairs sum n-14". Actual keys: {0,1} at 14 (sum 1), {0,2}{1} at 15 (sum 2). The sum is n-13, not n-14.
2. Claim (3) says 4056 "pairs" with sum n-13. Orbits were walked only at n = 13, 14, 15. Beyond 15 only the key was computed, which is a necessary condition for pairing. The claim should read "the key predicts, and orbits confirm to n = 15". The 4056 "onset" and the formula are stated as a law of orbits for n up to 30 and 20; they are a law of keys.
3. "Same phenomenon" as the H-020 failures rests on 2 of 6 words, with the other four not reached. That is honestly said, but the "now with a cause" sentence is stronger than the evidence. "Onset length" is a description of the key data, not a mechanism: nothing says why the onset for `4056` is 13 or for `45` is 8.
4. "13 cases compared, Coxeter partition = orbit partition" counts n = 13, 14, 15 rows of the same three cores and two others. That is about five independent cores and is thin support for "orbit = key class" as a working rule; the Next section proposes it as a prefilter, which is only safe in the non-pairing direction.

## New?

`grep` of F-053, H-021, H-020, E-052 in research/: F-053 already states that `3346` fails (five orbits at n = 14) and `4056` "pairs 0 with 1 and leaves 2 alone". E-051 (EXPERIMENTS.md l.167) records `350066` and `6600066` as the slide-of-one-or-two failures. Nothing found using the Coxeter key to partition offsets (only search gates and verifyMove checks, E-033 and neighbours). So new: the key-based proof that these are true non-equivalences rather than walk failures, the 4056 formula, and the onset reading. Not new: that the two are exceptions.

## Evidenced?

Mostly. The table gives sizes, closure, and both partitions, and the scripts are saved. Gaps: the 585-word scan output is described only in aggregate (309 / 7 / 123) with no stated criterion for "one pairing centre" beyond the sentence, and the cited long list is a script, not results. The key partition is stated for `3346` "to n = 40" and `3345` at four lengths only. The 585-word scan is at a single n = 24, so "onset" behaviour is untested across n for anything but the five named cores.

## Required for acceptance

1. Correct `6600066` to sum n-13 (or show n-14 at some length).
2. Separate orbit-verified statements (n <= 15) from key-only statements (beyond) in the claim, especially for 4056, and drop or soften "with a cause".
3. Either walk one 4056 orbit at n = 16 ({0,3},{1,2}) to confirm pairing there, or state that pairing beyond 15 is unconfirmed.
4. Save the scan's summary output (counts and the list of 7 triple-class words) in the round directory or the submission, so the 309/7/123 can be checked without a 3-minute run.
