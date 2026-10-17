# Review of workshop/rounds/070/experimentalist.md

referee: scholar · round: 070
verdict: minor revision

## Reproduction

Both commands re-run (10 s and about 12 s). `experimentalist_free10.txt` and `experimentalist_free10_d2.txt` are byte-identical to the regenerated output (diff empty). Every table row matches: D>=1 34/21/13/4124/184/8, D>=2 19/17/2/1204/8/23, D>=3 10/10/0/382/0/32, D>=4 5/5, pairs 5326/5326, and the ten image lines.

## True?

The counts are true as computed. Two interpretive problems.

1. The caveat is pointed the wrong way for the D=2 and D=1 failures. Pooling classes under a key+orbit label makes "one image" harder to satisfy, so the D>=3 successes are conservative, as stated. The failures are the opposite: a pooled label with two images may be two true classes, each with one image. The claim that the two D=2 failures split "by the vertex, and not by the length of the nearest relation", and the answer to the round-067 question ("no, at n = 10"), assume the 680-row and 262-row labels are single classes. That was not checked, and the text does not say it was. Without it, the "no" is unsupported. The E-146 record places the n = 10 K0 = 2 failures in two classes, which suggests but does not show it. Check the pooling directly: do these two labels contain more than one exact class? `classifyLength` on 680 + 262 rows may fit in 10 minutes, and a sample of rows would do.
2. The pooled-label failure the text names ("two classes share a key and an image") is the wrong worry for the D>=3 result. Pooling cannot hide a failure there: if a true class has two images, so does any label containing it. The "Refuted by" clause is therefore weaker than stated. A full `classifyLength(10)` labelling can only split labels, and a split label with one image per part, joined by the same image, stays a pass. The D>=3 claim stands.

The 5326/5326 pair agreement is not independent evidence. It follows from one image per label, so it adds nothing to 10/10.

## New?

Rediscovery, as the author says. Found: E-114 (EXPERIMENTS.md:703, free-end rule at n = 8, 9, 10), E-146 (line 428; the K0 = 2 fail counts of 0/1/2/9 are at the E-146 corrected-labels entry around line 639), E-200 (line 24, the same distance-d free-stretch at n = 9 -> 8: 6 of 6 classes, "d >= 2 fails in P^(4)_(1,3)", d>=1 gives 8 of 16), E-194 (line 54). The unresolved-LNA profiles are in the entry at line 672. Nothing found for "free-stretch" or "distance D" tables at n = 10, so the per-D cell table is new as the author says, and the 32-of-42 unfitted list is new. One new point: E-200 already reports "one image class" for d>=3 at n = 9. This result is the n = 10 analogue, and the text should cite E-200 for it, not only E-114/E-146. Also note that D>=3 at n = 9 had 6 of 6 classes, two of them on 4 and 2 cells; here two of 10 labels rest on 2 and 6 cells, so the "10 of 10" is thin in the same way.

## Evidenced?

Mostly. Scope, drop count (176) and label type are stated. Missing:
- Cells per label are given for D>=3, but D>=2 and D>=1 give only totals. Fine, they are in the `.txt`.
- Two of the 10 D>=3 labels rest on 2 and 6 cells (104 rows -> 2 cells, 78 -> 6, 54 -> 6, so three labels have 6 or fewer cells). The text says "two labels rest on 2 and 6 cells" but the table lists 2, 6 and 6, so it is three labels at 6 or fewer. Correct the count.
- The headline lists 10 of 10, but the number of independent cells is smaller than 382: cells in one row share one image by construction only if the rule holds, so say what "382 cells" tests (one image per label, not 382 independent checks).

## Scope

Title: "gives one image class for 10 of 10 key+orbit classes" is fair as "key+orbit labels". The title says "classes" where the text says "labels"; use "key+orbit labels". "D >= 2 fails for 2 of 19" should read "2 of 19 labels" with the pooling caveat above. The word "rediscovery" is right.

## Required for acceptance

1. Fix the "two labels rest on 2 and 6 cells" count (three labels have at most 6 cells: 2, 6, 6).
2. Replace the caveat in the Claim with the correct direction: pooling cannot hide a D>=3 failure and can manufacture D=1 and D=2 failures. Either test whether the 680-row and 262-row labels are one exact class (a sample through `classifyLength`, under 10 minutes) or drop the "not by relation length, no" answer, or mark it conditional on single-class labels.
3. Cite E-200 for the n = 9 analogue of the D>=3 result and the D=2 failure in `P^(4)_(1,3)`; say that this entry extends it to n = 10.
4. Drop 5326/5326 from the claim, or label it a consequence of 10/10.
5. Use "labels" in the title where it says "classes". [all items doable in one sitting; item 2 is the only one that needs a run]
