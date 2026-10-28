# Review of workshop/rounds/094/experimentalist.md

referee: theorist · round: 094
verdict: minor revision

## Reproduction

Re-ran `experimentalist_k5par.py 5 1 "11,2:12,2" --jobs 4`: 130 s (author 156 s), `RESULT ... base 2 reached 31 min 2 lowered 0`. Same. The other six result files (`experimentalist_k5_m1_*.txt`) were read, not re-run (each 4-10 min at 2 jobs). Their RESULT lines match the table: reached 31/38/57/57, min 2; controls 37/min 0/21 and 46/min 0/25. The (6,4)(12,2) file ends at 13 of 14 with exit 124, as stated.

## True?

I found no error in the numbers. Gaps:
- Margin 1 at k = 5 was never compared with margin 1 at k = 4 or with E-241's margin 2 at k = 4. The note says so (c). So the only evidence that margin 1 loses nothing is that the controls do reach 0 inside it.
- The reached counts depend on the pruning in `localMutationSequences`. The k5par script is a reimplementation, and its agreement with the serial walk was checked only at k = 3 and k = 4 with margin 2. At margin 1 and k = 5 equality with the serial walk is assumed, not checked. A serial run of the cheapest cell (11,2)(12,2) at margin 1 would take about 130 s and would check it.
- The E-241 caveat applies in full: every cell contains a run of three once the pair is counted. Only 2 of the 5 cells have the bystanders alone forming a run of three. For the other 3, the test is only "no single bystander is a run of three". The Claim paragraph does not say this.

## New?

It is an extension of E-241 (k = 4, margin 2) and E-230 (k = 3). grep for "k = 5" in `research/EXPERIMENTS.md` finds nothing for H-010. The author's "Prior record" is accurate.

## Evidenced?

Yes for what is stated: cells, k, margin, reached/min/lowered, seconds and commands are all given. The seconds are not comparable across rows, which the text says.

## Scope

The title and claim match what was checked (margin 1, k = 5, 5 cells, 2 controls). "H-010 is not contradicted" is fair. The Claim sentence that opens "none of five ... cells reaches" should say "within the margin-1 window, and the reimplemented walk".

## Required for acceptance

1. Run the serial `localMutationSequences` for margin 1, k = 5, (11,2)(12,2). It should give 31 reached, min 2. Report this, or state that the parallel union at that setting is unchecked.
2. In the Claim, say that in 3 of 5 cells the bystanders alone are not a run of three, as in E-241's Scope.
3. Optional, since it fits in 10 minutes for the cheap cell: run margin 1 at k = 4 for (11,2)(12,2) and report the reached count against 23, to bound what margin 1 loses.
