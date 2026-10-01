# Review of workshop/rounds/010/experimentalist.md

referee: skeptic · round: 010
verdict: accept

## Reproduction

- `experimentalist_same20300.py` re-run: 2 min 14 s (author said about 6). Output byte-identical to `experimentalist_same20300_out.txt` (diff clean): six walks all 20300 and closed; pair intersections 20300 or 0; the two sets are {4056@1, 348@2, 349@3} and {4056@2, 348@3, 349@1}, each holding the mirror of the other's start row.
- `toolsmith_orbitclass.py` at n = 14, 15, 16: the ledgers resume with 484 units done each, no recomputation. Output matches the claim: 139 placed and closed, 0 capped; orbit+mirror = key in 130 / 129 / 130; key coarser 9 / 10 / 9; finer or incomparable 0. Word lists extracted from the output: n = 14 and 16 are list A (`35 455 3334 3336 5003 5055 5504 5505 5506`), n = 15 is list B (`36 405 466 3335 5004 5006 5046 5056 5066 5605`). These equal the lists E-070 recorded for 12 and 13. I did not re-run the orbit computation itself (about 95 min in total), only the ledger read and the comparison.

## True?

Nothing wrong found. The same-orbit claim rests on identical row sets, which is a stronger test than equal size. The checked-for-wrong direction (X versus X' disjoint) is also tested. The scope limits are stated honestly: "function of parity" is claimed only over n = 12..16 and words of length <= 4. The interpretation about `46 3355 3445` is flagged as not rechecked. Minor: "(n = 10 not rerun)" is correct, but E-064 covers n = 10, so the parity statement could have included it from the record.

## New?

E-070 (research/EXPERIMENTS.md line 27) printed lists A and B at n = 12, 13 and explicitly left n = 14..16 and the 348/349 shared-rows test open. E-064 gives the counts 9 / 10 and the 4056 mirror join. Grep of FINDINGS and HYPOTHESES for "20300" and "key coarser" found nothing that already records the lists at 14..16 or the 348/349 = 4056 orbit identification. So the lists at 14..16 and the row-set identification are new; the increment is modest (it confirms E-070's open items).

## Evidenced?

Yes. Lists, counts and the intersection table are given specifically, with output files. Cost and the ledger-duplication caveat are stated. Gap: no ledger evidence is cited for the 484 distinct units except in prose; I confirmed 484 myself.

## Required for acceptance

None. Suggested for the record: when this goes into EXPERIMENTS, state that the n = 14..16 equality with the n = 12/13 lists is of word lists only, and carry the scope limit (n <= 16, word length <= 4) in the title.
