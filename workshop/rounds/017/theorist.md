# `3334` and `2455` are one double mutation from `35`, so they lie in the class of `333@1` (not `333@0`, the `444` orbit); the same move sends `444` to `34`

author: theorist · round: 017 · kind: result (lemma checked by exhaustion over a bounded range, not proved)
thread: T1/T2/T4 · bears on: H-020, H-021, F-032, F-051, E-065, E-079, E-083, E-086

## Claim

1. **Lemma R (hand-derived, checked).** Take four consecutive relation starts with lengths `a, b, b, d`
   (3 <= a <= b, b >= 3, d >= b; relations as intervals `(s-1,s-1+a), (s,s+b), (s+1,s+1+b), (s+2,s+2+d)`).
   The (dual) double mutation of arXiv:2310.08346 at the middle relation `r = (s,s+b)` gives `(a-1, b, d+1)`; with no fourth
   relation, `(a,b,b) -> (a-1,b)`. Reason: a relation ends at `t+1` (the third), none ends at `s+1`; the new relation
   `(s-1,t-1)` swallows the shortened first one, the third is lengthened to start at `s` and so contains `r` (dropped), the
   fourth starts properly inside `r` and is lengthened to start at `s+1`. Instances: `444 -> 34`, `455 -> 35`, `3334 -> 2 3 5 = 35`
   (the 4 becomes a 5 because `r` covers its start), `4445 -> 346`, `4444 -> 345`.
2. **Class label.** For a closed orbit write `J` = the offsets `o'` of `333@o'` it holds (these come in pairs `{j, n-6-j}`; `333` is
   self-dual). At n = 12..16 `J` is `{0,n-6}` for the `444` orbit and for every placement of `34`, `444`, `2234`, `2444`, `4445`;
   `J = {1, n-7}` (orbit sizes 447/320/763/516 at n = 13/14/15/16) for `35`, `455`, `3334`, `2455`, 55.
   So `3334` and `2455` are in class 1, never class 0: at odd n all placements, at even n the placements of one parity
   (the others are in the class-less parity shadow, as for `35` itself).
3. **Reading.** Class of `3x` is `x - 4` for x = 4..7 (folds by `j <-> n-6-j` at x = 8, 9), at every placement (odd n). The "4" of
   `3334` is irrelevant: R turns it into the 5 of `35`. The 4 of `444` has nothing after it and R gives `34` (class 0).
   So "the `444` orbit" is the class of `34`/`333@0`, the small orbit is the class of `35`/`333@1`; E-086's `235/255/455`, `2455`, `3334`
   having equal size is one orbit, not coincidence.
4. **`k(33x) = 2x`.** Given the class table (`33x@o` has `J = {x+o-3, n-3-x-o}`, drift label of E-065, and distinct `J` are
   distinct orbits), `33x@o` and `33x@o'` share an orbit iff `o' = o` or `o + o' = n - 2x`, so `s = n-2x`, `k = 2x`, `d = x-3`.
   This is E-065's argument; what I add is the *link*: the class of `3x@0` is reached by one anchored rule `3x@0 -> 33(x-1)@0`
   (seen at n = 13 for x = 5 on a labelled path; consistent with the class table for x = 4..7), so the `3x` family is
   `x - 4` and joins the `33y` chain at its drift label.

Not claimed: that class 0 and class 1 are inequivalent derived classes (closed under this move set only; no invariant is known
and the P/Q nulls stand); the upper bound in 4 (orbit contains nothing beyond the chain) is enumerated at n = 12..16, not proved;
`k(33x) = 2x` is not newly derived, only its end link is located. Words with letters >= 6 and parity shadows are not explained.

## Evidence

- Lemma R tested at n = 14 against `doubleMutation.rewritesOf` at every interior placement: four-run 250 of 250 (a in 3..b, b in 3..7,
  d in b..9); three-run 119 of 119 (a in 3..b, b in 3..8); n = 12 three-run 77 of 77. The cases `a = 2` give a one-arrow
  relation, not an LNA row, and are not tested (105 "missing", all of them `a = 2`).
- Every placement of `3334` has `35@(o+1)` as a one-step double-mutation neighbour, `2455@o -> 35@(o+1)`, `455@o -> 35@o`
  (n = 12, 13, 14, 16: all placements; see `theorist_label_n*.txt`).
- Orbit classes by `J`, n = 12..16, closed orbits, full tables in `theorist_label_n{12..16}.txt`. Excerpt (n = 13, odd, no parity split):

| word | class J | orbit size | placements |
|---|---|---|---|
| 34, 444, 2234, 2244, 2444, 4445 | {0,7} | 2386 | all |
| 35, 55, 235, 255, 455, 2455, 3334 | {1,6} | 447 | all |
| 36, 466 | {2,5} | 449 | even o (odd o: 224 shadow) |
| 37, 77, 237, 277, 2477, 3336 | {3,4} | 501 | all |

  At n = 14, 15, 16 sizes of class 0 / class 1 are 3767/320, 5648/763, 8134/516 (class 1 at n = 14, 16 is one parity of the
  placements; other parity is a 272 / 446 orbit with no `333`).
- Labelled shortest path (n = 13) `3334@3 -> 3@2`, 10 steps: two double mutations to `5@..` and `35`, then `35@0 -> 334@0` (anchored rule
  `3 5 -> 3 3 4`), `334@0 -> 333@1` (x = 4 rule), `333@1 -> 3@2` (`333 -> 23` at w = 2). `444@3 -> 3@1` and `-> 4@0`: 5 steps, `444 -> 34` first.
  `4@0 -> 3@1` is in class 0 (E-080).
- Hand check of R on `3334` placed at starts 4..7: intervals `(4,7),(5,8),(6,9),(7,11)`; `r = (5,8)`, third `(6,9)` ends at `t+1 = 9`, none ends at 6;
  result `(4,6),(5,8),(6,11)` = `2,3,5`; strip the 2: `35`.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/017/theorist_rrule.py 14     # lemma R, ~3 s (also 12)
for n in 12 13 14 15 16; do timeout 10m .venv/bin/python workshop/rounds/017/theorist_label.py $n; done   # classes, 1..25 s each
timeout 10m .venv/bin/python workshop/rounds/017/theorist_class.py 13 34 35 36 37 55 235 455 2455 3334 444 2244   # J of listed words
timeout 10m .venv/bin/python workshop/rounds/017/theorist_chain.py 13 3334 3 3 2     # labelled path
```

## Prior record

E-083/E-086 record the memberships (by size, then by row set) and call the cause open; E-065 gives the 33x drift and says
the end link is open; E-080 has `4@0 -> 3@1` as one width-4 move; E-074's list A contains `35 455 3334`. New: the one-step reduction
to `35` and the shared lemma R with `444 -> 34`, and the class label `J` as a table over words (not in the record; grep
`J`, "class" and `35` in `research/` found nothing equivalent). Not in `RETRACTIONS.md`.

## Code changed

None.

## Next

- experimentalist: run `theorist_label.py` over all 4-letter words with a 4 at n = 13, 15 (cheap) and check whether the `J`
  predicted by "apply R, then read `3x -> x-4`" matches for every merged word; any miss is a counterexample to item 3.
- skeptic: item 3's `x - 4` fold at x = 8, 9 and n = 17 `3334` (one orbit size, ~minutes) as a test; is the claim "class 1 != class 0"
  anything but a statement about this move set?
- theorist (next): write the anchored `3x@0 -> 33(x-1)@0` as a rule-table row and derive `J` for `3x` at all x from it and R.
