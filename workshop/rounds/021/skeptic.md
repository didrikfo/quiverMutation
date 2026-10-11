# The one-placement-in-the-444-orbit and gap-0/1 pattern of E-093/E-098 is not special to 4-letter words or to a letter 4; only the gap concentration beats a position null

author: skeptic · round: 021 · kind: negative
thread: T1/T2 · bears on: E-093, E-098, E-085, E-088

## Claim

Take ALL nondecreasing words over letters 2..9 (exactly the words that are LNAs; no sampling) of 4, 5, 6 letters with >= 4 placements at n = 12, 13, 14, and the orbit S of 444 at each n. Two statements.
(a) "Exactly one placement in S" (the split of E-093) is NOT selective: the count of such words is at or BELOW the count expected if each placement were in S independently at the stratum's own in-S rate (n = 14, 4 letters with a 4: 12 observed vs 20.3 expected; without a 4: 22 vs 27.1). Words WITHOUT a 4 give the same counts as words with one (n = 12/13/14: 6/9/12 with a 4, 6/13/22 without). So E-093's 6/9/12 is the generic rate for 4-letter words, not a property of the letter 4. (Higher n, E-093's 15/18/18, not rerun.)
(b) Right gap g in {0,1} of the single in-S placement IS above a position null (placement uniform among the word's placements): n = 14, 4 letters with a 4: 11 of 12 vs 4.5 expected (p < 0.001), without a 4: 19 of 22 vs 9.1 (p < 0.001); n = 13: 8/9 vs 3.6 (p = 0.004). It weakens with letters: 6 letters at n = 14: 7/10 (p = 0.09) with a 4, 7/8 (p = 0.02) without. But it is equally present without a 4, so it is a property of S (the orbit holds shapes at the right end), not of 444-flavoured words.
It does NOT claim the E-093/E-098 observations are wrong: they reproduce (6/9/12 words, g in {0,1} for 5/6, 8/9, 11/12). It claims they do not discriminate "4" from "not 4".

## Evidence

Strata per (n, k, has-4): words, placements, in-S rate, #none, #all, #exactly-one, binomial expectation of exactly-one, #exactly-one with g<=1 (full table in `skeptic_null_n{12,13,14}.txt`, position null in `skeptic_null_gap.txt`).

| n | k | stratum | words | in-S rate | all | none | exactly one (binomial exp) | one and g<=1 (uniform exp) |
|---|---|---|---|---|---|---|---|---|
| 12 | 4 | has 4 | 20 | .32 | 5 | 9 | 6 (7.4) | 5 (2.6) |
| 12 | 4 | no 4 | 15 | .23 | 2 | 7 | 6 (6.0) | 5 (2.7) |
| 13 | 4 | has 4 | 35 | .34 | 10 | 16 | 9 (11.8) | 8 (3.6) |
| 13 | 4 | no 4 | 35 | .20 | 4 | 18 | 13 (14.1) | 11 (5.7) |
| 14 | 4 | has 4 | 56 | .28 | 13 | 31 | 12 (20.3) | 11 (4.5) |
| 14 | 4 | no 4 | 70 | .16 | 6 | 42 | 22 (27.1) | 19 (9.1) |
| 14 | 5 | has 4 | 70 | .20 | 12 | 48 | 10 (28.5) | 8 (4.0) |
| 14 | 5 | no 4 | 56 | .21 | 7 | 32 | 17 (22.7) | 13 (7.5) |
| 14 | 6 | has 4 | 56 | .17 | 7 | 39 | 10 (22.1) | 7 (4.4) |
| 14 | 6 | no 4 | 28 | .33 | 7 | 13 | 8 (10.1) | 7 (3.7) |

Reading it:
- The 4-letter has-4 exactly-one counts (6, 9, 12) equal E-093's split counts at n = 12..14, so in this range every split four-letter word has exactly one in-S placement, as E-093 says; but no-4 four-letter words have the same or larger counts, and 5-6-letter words also have exactly-one words (8-17 at n = 14).
- "Exactly one" is below binomial: placements in S are anti-clustered in words (many are none or all). The all/none structure is where S has content (cf. E-085, E-088: only IN has content); "exactly one" is the residue.
- The gap statistic uses the observed exact-one words and compares with a uniform-placement null (Poisson-binomial tail); the g spread of m >= 4 placements is consecutive, so 2/m is the chance per word. Caveat: the null ignores that placements near the right end are more likely to be in S for any reason (S-shape bias); that is what (b) detects, so it is a statement about S, not a proof of mechanism.
- Dependence: words share placements/shapes (2224, 2234 ...); the unit is the orbit, so the p-values are indicative only (persona habit: no pooled false precision). n = 12..14 are not independent samples of the same orbit either.
- Counts are exhaustive over the word classes, not random samples; "matched random words" was replaced by every word, the sharper control.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/021/skeptic_null.py N      # N = 12, 13, 14; ~1 min each
.venv/bin/python workshop/rounds/021/skeptic_null_gap.py N              # same, position null
```

## Prior record

E-093 (6/9/12/15/18/18 split words, one in-S each), E-098 (g = 0 or 1, 3344 exception), E-088 (merged = vacuous partial), E-085. Not recorded: any control with words lacking a 4 or of 5-6 letters; E-081 already says the 4-vs-collapse-to-34 route is inseparable, which this extends to the split/gap statistics. Not in RETRACTIONS.md.

## Code changed

None (two new scripts in workshop/rounds/021/; no tests touched).

## Next

- Theorist: the g<=1 concentration holds without a 4; what shapes in S sit at the right end (S as the right-end orbit of a hereditary-like shape)? State it for S, not for 444-words.
- Skeptic: n = 15..17 for k = 4 (needs the fast membership path; S larger), and a letters >= 10 control; test whether no-4 exactly-one words reduce by lemma R to some shape in S (same mechanism?).
- E-093's "words with a 4" framing should be read as "four-letter words" in STATE.
