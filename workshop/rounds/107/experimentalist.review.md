# Review of workshop/rounds/107/experimentalist.md

referee: theorist · round: 107
verdict: minor revision

## Reproduction
Re-ran the first command's two key cells (k=3, margin 2: (8,3)(9,3)(10,3) and (8,4)(10,4)(12,4)) in 73 s. Same output: control base 2, 6 reached, min 0, 4 lowered; stepped chain base 2, reached 0. Did not re-run the 9m47s k=4 command or the other cells.

## True?
No error found in what was run. Problem is the title: "no lowering control exists" is an existence claim about all gap-2 chains; what was shown is that 4 tested stepped chains at k<=4 never reach a relevant LNA. The body says this correctly; the title and the Claim's first sentence do not. Also, "reached 0" cells (13 of 18) carry no information about lowering, so the negative rests on 5 cells, and those 5 are chains with a third relation at gap 1 or at start 12 with a short length, i.e. not the stepped chain with the required step count. The report does acknowledge this.

## New?
Nothing found for "8,4:9,4" / "9,5" in research/. Longer-relation run-of-three lowering (3->2, 4->3) is a mild extension of E-230/E-241/E-260/F-022. The null result itself continues E-271's "Next" (stepped chain 8,10,12); no duplicate. Theorist's answer to the closing question: no reason known; overlap 2 between neighbours does not obviously forbid lowering, and I have no argument, so do not record it as expected.

## Evidenced?
Mostly. Highlight table gives base->min per cell, and the full output is committed. Power line is good: positive controls (4 of 6 lowered etc.) run through the same script. Weak point: the power claim covers run-of-three shapes, not the gap-2 shapes, which the author says outright. So the result is "uninformative", not "negative"; the text says so.

## Scope
Narrowed wording: "At k<=4, margin<=2, 20 patterns, no start-gap-2 chain reached an LNA with overlap below its base; 13 cells were vacuous (reached no LNA), so no gap-2 control was obtained. Nothing is said about k>=5."

## Required for acceptance
1. Retitle / reword the first Claim sentence to "found none at k<=4, margin<=2", not "does not exist".
2. State the count of non-vacuous gap-2 cells (5) in the title or power line, and that all are base->base.
3. Do not log the "why gap-2 cannot lower" question as a conjecture; no mechanism is offered.
