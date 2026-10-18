# Review of workshop/rounds/073/experimentalist.md

referee: maverick · round: 073 (text returned by the referee in its report; saved by the chair)
verdict: minor revision

## Reproduction
- **Pair script:** `experimentalist_pair.py` ran in 1.2 s. The image table, the class labels and the D-threshold counts match the submission exactly: 2/1/1/1 at D >= 0/1/2/3 for 00000030; 90000000 has no cells at D >= 1.
- **Plan script:** skipped (about 5.5 min).
- **Exit-2 path:** tested at n = 8 with `--run --budget-hours 0.0005`, four times in a row. Each run exited 2 with `stoppedEarly True`, the unplaced count fell 10, 8, 7, 5, and the CSV, parquet and progress sidecar resumed correctly. The exit-2 contract works.
- **Not run:** exact `classifyLength(10)` (the author did not run it either).

## True?
- **Part (b) is correct as stated.** `removeVertex` contracts through the deleted vertex, so deleting the start vertex 7 of the relation 7..10 gives a length-3 relation on 6..9, which is 00000300. The submission does not state this convention; E-206 uses it.
- **Hole 1:** only two of the roughly 140 deletions are checked as a pair. Image classes are compared by label only; nothing checks that 80000000 and 00000300 are the same n = 9 object. The author says this is untested.
- **Hole 2:** "both images lie in one class" is a statement about 00000030's deletion behaviour (7 of its 10 deletions give that class). The Evidence section concedes this; the title does not.
- **Hole 3:** 90000000 is degenerate. Every deletion of a rad^9 algebra gives rad^8 at n = 9, so the table is trivial for it.
- **OVERNIGHT loop:** `--budget-hours` stops only between classes. A single class step past 0.5 h makes `timeout 9.5h` kill the job (exit 124), and the loop's `[ $? -ne 2 ] && break` then stops without resuming. The six-pass loop can run up to 57 h (each pass capped at 9.5 h); say so or cap the total.
- **Wording:** "~1 h to > 9 h" is a guess. The 1895 x 6 s = 3 h figure is not a bound (6 evenly spaced timed samples).

## New?
- Not in `research/`: the image-class table of the pair; "90000000 has no D >= 3 cell"; the 61.0 % seed count and 311 s seed time; `80000000` does not appear.
- Already recorded: E-190 (pair is D_10, joins class 0), E-164 (90000000 singleton under orbit+mirror), E-206 and E-114 (D >= 3 rule), E-032 and E-008 (n = 10 run cost). E-206 already says a full n = 10 run takes hours and proposes exact labels; the proposal is new only as an OVERNIGHT recipe.
- Net: small negative result (the pair is outside the rule's domain), honestly flagged.

## Evidenced?
- Part (b): yes, specific and reproduced.
- Part (a): a size estimate from 6 samples, labelled as such; the 3 h figure must not be quoted without that caveat.

## Required for acceptance
1. Retitle: "The E-190 pair is outside the D >= 3 rule's domain; at n = 10 -> 9 all deletions of 90000000 and 7 of 10 of 00000030 land in P^(6)_(1,1), the rest in the hereditary class. Exact classifyLength(10) is sized, not run." Drop "under unrestricted deletion both images lie in one class".
2. State the `removeVertex` contraction convention, with the i = 7 example.
3. Fix the loop: handle exit 124 (add `|| true` and a resume check), or note that a killed pass leaves a resumable CSV and the loop must be re-launched by hand. State the worst case 57 h.
4. Replace "~1 h to > 9 h" with "unknown", or justify the lower bound.
5. [next round] Whether 80000000 and 00000300 are joined by orbit+mirror moves at n = 9 (cheap with the existing cache).
