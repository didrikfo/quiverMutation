# Review of workshop/rounds/018/skeptic.md

referee: theorist · round: 018
verdict: minor revision

## Reproduction

Re-ran `skeptic_rowset16.py N 4 0` for N = 14 and 16 (seconds each). Output files byte-identical to the saved `_fast.txt` (diff empty). Recounted all six `_fast` files with grep: IN/OUT/PARTIAL = 5/9/6, 10/16/9, 13/31/12, 16/53/15, 19/83/18, 19/83/18 for n = 12..17, matching the table. Every PARTIAL line has inS = 1 at every n. `3334`, `2455`: OUT at every n = 12..17. Recomputed the word-list size independently: 165 words with a 4, of which 120 have >= 4 placements at n = 16, 17 and 20, so the plateau at 120 is real. Walked file (`skeptic_rowset16_n16.txt`): 3334 and 2455 both show orbits (446, 516), overlap 0, closed. Did not re-run the 10-minute walk.

## True?

The headline holds as computed. Small errors:
- The walked run has 55 word lines (56 with the header), and the tags there are 13 IN / 28 OUT / 14 PARTIAL. "First 56 words, up to 3466" is loose by one. The walked-versus-fast agreement is for 55 words, about 46% of the list, not all of it.
- "offset is last or second-to-last for 17 of 18": checked against `skeptic_partial_where_n16.txt` and it holds (`3344` at 1 is the exception). This is a description, not a rule, and the report says so.
- The "list stops growing at 120 since letters <= 9 and the count of placements >= 4 bounds the sum" is a hand-wave. The recount above shows 120 is the number of words with >= 4 placements once n >= 16. The reason is the word's own length, not a bound on the sum. Say that instead.
- The remark that E-088 was vacuous for merged words is correct: orbits are disjoint, so "all IN or none" is a tautology for a word whose placements share one orbit.
- The 446-orbit is identified with E-090's 235-type orbit by size only. The report says so.

No counterexample found to any stated number.

## New?

E-088 (research/EXPERIMENTS.md:45) states "0 partial" for merged words only, so the all-words reading is not recorded there. E-085 gives the in-S counts, and the 5/10/13/16 consistency check is correct. E-090 (line 27) already explains `3334` and `2455` as `333@1`-class (lemma R), so their absence from S is explained there; the report cites it but presents it as a data point. grep for "partial", "split", "one placement" in FINDINGS/HYPOTHESES/RETRACTIONS found nothing on "exactly one placement of a split word is in the 444 orbit". That observation is new. Nothing in RETRACTIONS needs editing, but E-088's title wording ("0 partial") should get an erratum note, as the report itself suggests.

## Evidenced?

Adequately for the counts: ranges (n = 12..17, `--max-word 4`, letters 1..9, >= 4 placements) and the word lists are given. Missing:
- Whether "exactly one" is a theorem-shaped statement or an accident of the range. The data are 18 words at one n, with the same 18-word list recurring at 16 and 17 (12/15/18 earlier), so these are one family rather than independent cases.
- Words with 3 placements or fewer, words with letters >= 10, and words without a 4 are outside the range; "all words" in the title and claim should read "4-letter words with a 4 and >= 4 placements, letters <= 9".
- The other placements of split words were not classified (the report defers this).

## Required for acceptance

1. Restrict the claim's quantifier in title and Claim to words with letters <= 9 and >= 4 placements.
2. Correct "first 56 words" to 55 and state the fraction of the list walked.
3. Replace the "bounds the sum" explanation of the 120 plateau with the actual reason, or drop it.
4. State that the 17/18 "last or second-to-last" pattern is descriptive only. Better, give the offset against the number of placements for all 18.
5. Note that the 18 split words at n = 16 and 17 are the same list, so n = 16 and 17 are not independent evidence.
