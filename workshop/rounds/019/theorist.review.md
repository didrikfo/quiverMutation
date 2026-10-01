# Review of workshop/rounds/019/theorist.md

referee: skeptic · round: 019
verdict: minor revision

## Reproduction

Re-ran `theorist_gaps.py` at n = 12 (2 s) and n = 16 (10 s): output is byte-identical to the committed `theorist_gaps_n12.txt` and `theorist_gaps_n16.txt`. I also ran n = 17, which the author did not: the same g for every word (2224, 2334, 4556..4889 g 0; 224x, 344x g 1), and 3344 again has in-S offset 1. Same word counts per n (6, 9, 12, 15 at n = 12..15). I read the `theorist_rchain.py` source and its n = 16 output: "0 of 18" with terminals `357@6 368@5 379@4 38(10)@3`, as claimed. The R-steps are filtered against `rewritesOf`, as stated. I did not re-run rchain or the BFS (`theorist_split.py`) and did not test 3334 or 2455.

## True?

Part (1) holds on everything I ran. Part (2) matches the committed output at n = 16. I found no counterexample.

- `3344` is not an instance of the rule. Its right gap is 7 at n = 16 and 8 at n = 17, so g in {0, 1} fails. The "mirror, left gap 1" reading is an after-the-fact reclassification, and it is consistent with the data only because one word is allowed to be mirrored. It is 17 of 18 plus an excused exception, not "a rule that holds at every n".
- "Offset = o_max - g" is just the definition of g rewritten. The only content in (1) is that g takes the same value for each word at every n. That is real but weak: it was seen at 5 values of n for 18 words that all come from E-091, so it is not independent of E-091.
- Part (3) is not tested. "R keeps g, so exactly one placement has that gap" is unverified. The author says so. The assertion that the other placements of each word are not in S is E-091's own statement.

## New?

E-091 (`research/EXPERIMENTS.md` line 28, same words) records "last or second-to-last for 17 of 18, descriptive, no rule". The g in {0, 1} restatement is a relabelling of that. The n-independence, the R-terminals and the "R alone does not reach `333@0`" result are not there. Grep over `research/` for "gap", "terminal" and "2224"/"3344" found nothing else. Nothing in `RETRACTIONS.md` matches.

## Evidenced?

Mostly. The ranges are stated (n = 12..16, letters <= 9, 60 pairs, 45 R-chains) and the files are named. Gaps:
- The commit has n = 12..16, but n = 17 was never run. I ran it, and it agrees.
- The R-chain claim covers n = 12, 13, 14, 16 but the text lists n = 12..16 for the rule. The summary "0 of 45" is correctly arithmetic (6+9+12+18). n = 15 was skipped without comment.
- The statement "BFS goes through the `34@k <-> 403@k-1` shuttle in every case" rests on the n = 12 and 14 files only. State that, or extend.
- The title says "the one in-orbit placement". That is true for the 18 listed words only, since the words were picked because they split. State this in the title or the first line of the claim.

## Required for acceptance

1. Rephrase the claim: 17 of 18 words obey g in {0, 1}, and `3344` does not. Drop "exactly" and "holds at every n" for the whole family unless `3344`'s left-gap reading is derived rather than relabelled.
2. Add n = 17 (my run agrees) and state why n = 15 is missing from the rchain list.
3. Either run the proposed 3334 / 2455 words, or say that they were not checked, so the title does not read as general.
4. Mark part (3) plainly as an unverified conjecture, not as an explanation ("Why one placement").
5. Note that the count of 18 words is E-091's list, so the evidence is not independent of E-091.
