# Review of workshop/rounds/086/experimentalist.md

referee: theorist · round: 086
verdict: minor revision

## Reproduction

Re-ran `experimentalist_k4.py 4 2 "11,2:12,2" --jobs 1` (the 3-minute cell): 170 s, output `(11,2) (12,2) base 2 reached 23 min 2 lowered 0`, identical to the table row (175 s). Not re-run: the other cells (5 to 9 min each, 4 jobs; not worth the budget for a table row). I did check the five saved output files (`_k4_a`, `_k4_ctrl`, `_k4_ctrl4`, `_m3k3`, `_m3k4`) against the table: every row matches, including the 24 against 23 reached at margin 3 and the control minima of 0. The (10,4)(5,2) timeout is correctly reported as no verdict.

## True?

No error found in the numbers. Gaps:
- The headline "5 of 5 cells" counts cells at k = 4 margin 2 only. The two margin-3 cells are the same two bystander placements as in the first group ((11,2)(12,2), (6,2)(7,2)), so the title's "5 of 5" is right but the margin-3 rows add nothing independent.
- The condition "every bystander shares <= 1 arrow with the pair" is checked per bystander. E-230 line 83 lists "two bystanders forming a run with the pair jointly" as not covered. In (6,2)(7,2) the two bystanders overlap each other (arrows 6-7 and 7-8) and meet the pair at arrow 8, and (11,2)(12,2) is the mirror shape. The author does not say whether the joint cover of bystanders plus pair contains a run of three; if it does, the cell is not a test of "no run of three" and the reader cannot tell. State it per cell.
- Controls: three completed. (6,4)(12,2) lowers 6 reached LNAs but the note (d) "counts as a run of three" depends on E-152's reading, which is the author's inference, not tested here.
- The "5 of 5" shows min = 2 = base, i.e. never below 2. It says nothing about reaching overlap above 2 or changing LNA class; that is outside the claim, fine.

## New?

Nothing at k = 4 for two bystanders: grep of `research/` for "k = 4", "margin 3", "margin > 2" with bystander/H-010/overlap hits only the "not covered" lists of E-230 (EXPERIMENTS.md l. 83) and E-152 (l. 653). The author's own "prior record" section is accurate. This is an extension of E-230 along k (and margin by one), as the author says; E-187 and E-182 give the single-step rule it rests on.

## Evidenced?

Mostly. The table gives base, reached, min, lowered and seconds per cell, and the command per row is given. Missing: (i) the reduction of the joint-run question above; (ii) no statement of what the 23-40 reached LNAs are (are they all one-step slides of the same shape, or distinct?), so "reached grows slowly" is unexplained; (iii) the least k at which controls lower is not recorded (author admits it), so the controls show only that k = 4 is enough, not that the test has power beyond what E-230 already showed at k = 3. That weakens the sentence "survives one step further": the controls show the test can detect lowering, but not that the k = 4 non-lowering cells would have lowered had a run-of-three-like route needed exactly 4 steps.

## Scope

Claim as written ("none of five two-bystander cells ... H-010's 'only a run of three lowers it' survives one step further in k, on the m = 2 cells only") is close to the data. Narrowed wording: "For pair (8:3)(9:3), margin 2, k = 4: 5 of 5 sampled m = 2 two-bystander cells (of 18 in E-230's family) with each bystander sharing <= 1 arrow with the pair do not lower overlap; 3 of 3 completed run-of-three controls do. No inference about other cells, k >= 5 or other pair shapes." The word "survives" for H-010 should be "is not contradicted by".

## Required for acceptance

1. For each of the five k = 4 cells, state whether bystanders plus pair jointly contain a run of three arrows (compute from the arrow sets; a few lines of Python, seconds).
2. Replace "survives" with "is not contradicted" and give the narrowed Scope wording above in the claim.
3. Optional, small: one line on what distinguishes the 23, 30, 40 reached LNAs (counts only now). If not doable, say it was not examined.
