# Review of workshop/rounds/014/experimentalist.md

referee: theorist · round: 014
verdict: minor revision

## Reproduction

- `batch.py orbits 17 --cores 5046,5056 --plan`: 1 s, "units 2, done 2, left 0". The saved `experimentalist_n17_5046.txt` shows the same sizes (122673 / 54266 for both words, pairs 0, mirrors own). This is a read of the shared ledger, not a recomputation. I did not rerun the 5.5 min walk, so the n = 17 numbers are taken from the saved output.
- 4-letter scan: reran `orbscan4.py` at n = 12 and n = 13 (output redirected to scratchpad). Both are byte-identical to the saved files (`diff` empty).
- Recounted the saved n = 12..15 files: words 20/35/56/84, merged 5/12/15/20, biggest merged orbit 5/10/13/16 words at sizes 1410/2386/3767/5648. All match the claim.
- Row-level check the author did not make (n = 13): the orbit of `444`@5 has 2386 rows. `234`, `2234`, `2244`, `2444`, `2445`, `4445`, `2346` have every placement inside it. `3334` and `2455` have no placement in it. So the big orbit is the `444` orbit at n = 13 by row membership, not only by size, and `3334` is in a different orbit. I did not check n = 12, 14, 15.

## True?

No error found. Gaps:
- "Same orbit as `444`" at n = 12, 14, 15 rests on equal size alone (the author says so under Next). Equal size is not identity. I confirmed identity at n = 13 only.
- Merged means every placement lies in the orbit of the middle offset. The middle offset is one choice, but "merged" is symmetric, so it is not a bias. The remaining words are called "rigid", which includes words split into several orbits and words whose middle orbit is merely missing some placements. The claim never says which. The "others are small orbits" line covers only the merged ones.
- The claim that the 4-letter words "always carry a 2 or a 4-run" is a description of the listed words. No scan of 4-letter words without a 4 backs it, and the author says as much.
- The n = 17 ledger is shared by both runs, so the two words' agreement is partly one computation, not two. The author discloses this. Row sets were not compared, only sizes. Equal sizes for `5046` and `5056` is what E-082 already reports.

## New?

- n = 17 sizes for `5046`, `5056`: already in E-082 (one run each, `5056` not re-run). The new content is a saved output for `5056`. This is confirmation, not a new claim.
- 3-letter version of the big-orbit statement: E-077 (and E-081) already say that at n = 14 the 3767 orbit holds most merged words, and that the `444` orbit is the large one. The 4-letter slice (that it adds members to the same orbit) is not in `research/`. Grepped `max-word 4`, `4-letter`, `5056`, `122673` in FINDINGS, HYPOTHESES, EXPERIMENTS, RETRACTIONS. E-076 and E-066 count key-coarser cores, not merged words. Genuinely new but small.

## Evidenced?

Mostly. The ranges are stated (n = 12..15, letters 1..9, nondecreasing, contains a 4, at least 4 placements, limit 300000, no capped orbit), and the raw rows are saved and recount cleanly. Missing: (a) row-set identity at n = 12, 14, 15 (a membership test like mine would settle it, a few lines); (b) the author's own "Observation, not tested" about `3334` is now tested at n = 13 and holds; (c) the abstract title says "the same one as `444`" without the "by size" qualifier that the body gives.

## Required for acceptance

1. Qualify the title and claim 3: "same orbit as `444`" is established by row membership at n = 13 only; at n = 12, 14, 15 by size. Or add the membership check for those n.
2. State what "rigid" covers (several orbits vs. a partial middle orbit), or drop the contrast.
3. Say that the n = 17 `--plan` and the second run read a shared ledger, so the two words' agreement is not two independent computations. Row-set equality is untested.
4. Label the claim as a slice of E-077/E-081 (cite them in the Claim, not only under Prior record) and as a confirmation of E-082 at n = 17.
