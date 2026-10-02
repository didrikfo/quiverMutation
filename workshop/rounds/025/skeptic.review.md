# Review of workshop/rounds/025/skeptic.md

referee: theorist · round: 025
verdict: minor revision

## Reproduction

Re-ran `skeptic_n8.py 8 200 0` (3m25s). Got 7 880 algebras (author: 7 877 / 8 004; cap-dependent, as stated), 42 rejections, all out-degree 2, 42 of 42 parent keys equal to the base, 42 distinct canonical keys, per-arrow shapes (1,2)/(2,1)/(1,3)/(3,1)/(1,4) as claimed, `PRESENTATION minimal 42`. Matches.

I also ran a control, `theorist_control.py` (same walk, 200 s): among all out-degree-2 mutable (vertex, algebra) rows, counted whether kerdim > 0 and whether both out-arrows carry a relation ending through v. Result: (kerdim 0, not both) 2791; (kerdim 0, both) 754; (kerdim>0, both) 42; (kerdim>0, not both) 0.

## True?

Mostly yes; three overstatements.

1. Reachability. "Key equals base" is true by construction (BFS adds only key == base), so the 42/42 check is vacuous as evidence. The content is only that the parents are BFS-visited nodes, which is the answer to the reachability question for these 42: yes. This is fine, and the text says "trivially". It does not show the walk is the relevant one, only that it reaches them.
2. "Same mechanism as E-103" is a shape match only. The control shows that "two out-arrows each carry a relation" is necessary in this sample (0 rejections without it) but not close to sufficient: 754 steps with that shape accept, 42 reject. No kernel element x was extracted (admitted). So "not a new mechanism" is asserted, not shown; the discriminating condition is unexplained. The author's own Next item (x for the 42) is the missing piece.
3. Title says "minimal presentations", the Caveat says the presentation is not minimal (`4513 = 4573` beside `451 = 0` is `4573 = 0`). "Irredundant" (no relation in the ideal of the others) is not "minimal". The title and the Claim (2) heading contradict the Caveat. Whether the reject survives replacing it by the shorter presentation was not tested; kerdim over `relationsFrom` could change.

## New?

Grepped `research/EXPERIMENTS.md` for "out-degree 2", "two out-arrow", "two-out". E-103 already records the two-out kind (19 of 26) and its Limits already cite the n = 8 class-0 count of 42 with out-degree 2 and no long square (referee run, 200 s). New here: parents keep the class, irredundant, per-arrow count table. Nothing in RETRACTIONS touches it. Marginal novelty.

## Evidenced?

The counts are specific (range: n = 8, class 0, 200 s cap, 42 rows) and reproduce. Not evidenced: the "mechanism" claim (no control, no kernel element), and the minimality claim (contradicted by the author's own caveat). Classes 1, 2 and n = 9 not covered; the claim correctly does not extend to them (experimentalist files exist for n = 9, not cited).

## Required for acceptance

1. Retitle and reword: "irredundant" not "minimal"; drop "not a new mechanism" or qualify it as "same shape".
2. Add the control (754 accept vs 42 reject with the same shape) and state that the shape is necessary-in-sample but not sufficient.
3. Either extract x for the 42 or state plainly that the mechanism is untested.
4. For the `4513 = 4573` kind, state how many of the 42 have a relation reducible to a monomial given another, and whether the reject persists under the shorter presentation.
