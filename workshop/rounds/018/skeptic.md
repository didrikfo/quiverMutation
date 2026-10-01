# At n = 16 (and 12..17) 18 of 120 four-letter words with a 4 split across the 444 orbit; 3334 and 2455 stay outside it

author: skeptic · round: 018 · kind: result
thread: T1/T2 · bears on: H-021, E-086, E-088

## Claim

Let S be the row set of the closed orbit of `444` (|S| = 1410, 2386, 3767, 5648, 8134, 11340 at n = 12..17). Over all 4-letter nondecreasing words with letters 1..9 containing a 4 and at least 4 placements (the `--max-word 4` list: 20, 35, 56, 84, 120, 120 words), the placements of a word are all in S, all outside, or split. Split words exist: 6, 9, 12, 15, 18, 18 at n = 12..17 (n = 16: 2224 2245 2246 2247 2248 2249 2334 3344 3444 3445 3446 3447 3448 3449 4556 4667 4778 4889). So "every word with a 4 has all placements in S or none" is **false** as stated for all words; E-086 does not say it (it is about *merged* words). Every split word has **exactly one** placement in S (n = 14: 12 of 12; n = 16: 18 of 18), the offset is the last or second-to-last placement for 17 of 18 words at n = 16 (4556 at 6 of 0..6, 2224 at 8 of 0..8, 3444 at 7 of 0..8) and offset 1 for `3344`. `3334` and `2455` have no placement in S at n = 12..17 (all 5..10 placements OUT); at n = 16 their placements lie in orbits of sizes 446 and 516 (walked, closed, disjoint from S).

Not claimed: that the one placement in S is explained, or that the "all IN" list (19 words at n = 16: 2234 2244 2346 2444 2445 2446-2449 2456 2467 2478 2489 4445 4456 4467 4478 4489 4566) has a rule.

## Evidence

Referee point on E-086 (and my own r015 line "0 partial"): a *merged* word has all placements in one orbit, and orbits are disjoint, so "all in S or none" is true for merged words by definition. That check was vacuous; its only content is the IN-count. The nontrivial question is the one above, over all words, and it has a negative answer: a word whose placements are in different orbits (rigid) can have one of them be the 444 orbit.

Counts by n (fast membership run, `LIMIT 0`, no orbit walks; 12 s at n = 16):

| n | S | words | all IN | all OUT | split |
|---|---|---|---|---|---|
| 12 | 1410 | 20 | 5 | 9 | 6 |
| 13 | 2386 | 35 | 10 | 16 | 9 |
| 14 | 3767 | 56 | 13 | 31 | 12 |
| 15 | 5648 | 84 | 16 | 53 | 15 |
| 16 | 8134 | 120 | 19 | 83 | 18 |
| 17 | 11340 | 120 | 19 | 83 | 18 |

All-IN at n = 12..15 is 5/10/13/16, same as E-083/E-086's in-S 4-letter counts (5/10/13/16); this is a consistency check of the identity of S. Words with 4 or more placements only; the list stops growing at 120 since letters <= 9 and the count of placements >= 4 bounds the sum. The split words' single in-S offset is listed in `skeptic_partial_where_n16.txt`: for 4556, 4667, 4778, 4889 it is the last placement; for 2245..2249 and 3445..3449 it is the second-to-last. No rule is derived.

Walked (limit 300000) at n = 16: the first 56 words, up to `3466`, give identical tags to the fast run (the walked run timed out at 10 min on big outside orbits; the fast run needs none). Outside orbits of `3334`, `2455` at n = 16: (446, 516), overlap with S 0, all closed; these are the `J = {1, n-7}` orbits of E-088 (sizes 446/516 against E-088's 320/516/763 list at n = 14/15/16: 516 at n = 16 agrees; 446 is the 235-type orbit by size only, not compared as a set).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/018/skeptic_rowset16.py N 4 0     # N = 12..17, about 12 s each; output skeptic_rowset16_n{N}_fast.txt
timeout 10m .venv/bin/python workshop/rounds/018/skeptic_rowset16.py 16         # walked, stops at 10 min at word ~3466; skeptic_rowset16_n16.txt
.venv/bin/python workshop/rounds/018/skeptic_partial_where.py 16                 # skeptic_partial_where_n16.txt
```

## Prior record

E-086 (n = 12..15 merged words: 0 partial) and E-083 (in-S 4-letter counts), E-088 (`3334`, `2455` in class `J = {1, n-7}`, n = 12..17). New: the all-words reading, 6/9/12/15/18/18 split words, each with one in-S placement; the vacuity of "0 partial" for merged words; n = 16, 17 counts. The n = 17 row (S size 11340) agrees with E-088's referee value for the `444` orbit.

## Code changed

None in the library. New scripts `skeptic_rowset16.py`, `skeptic_partial_where.py` in `workshop/rounds/018/`. No tests touched.

## Next

- theorist: why exactly one placement of a rigid word lies in the 444 orbit (e.g. does `4556@6` reduce by lemma R / the `3x@0` link to `333@0`?)
- skeptic (next): which orbit the other placements of split words occupy and whether the split words' remaining placements form one orbit (a "two-orbit" word) or many; and whether the all-IN list is exactly "words with a collapse path to `34`".
- E-086's wording should say "merged" tests carry no row-set content; the content is the IN count.
