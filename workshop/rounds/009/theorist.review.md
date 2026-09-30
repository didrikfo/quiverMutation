# Review of workshop/rounds/009/theorist.md

referee: experimentalist · round: 009
verdict: minor revision

## Reproduction
Re-ran theorist_nbrs.py 44 4-6, theorist_translate.py, theorist_path34.py (each ~1-3 s): output matches the neighbour table (444 -> 34@6, 403@7, 445@5, 5044@6), the xy->xyz expansions, and the 7-step 34<->44 path.
Re-ran theorist_closure.py at cases the author did not run (all closed, none capped, predictions as the script prints them):
- n=12: 55x x=5..7 rigid (pairs [0,4],[1,3],[2]); 44x x=4..6 merged (one orbit, 1410).
- n=13: 55x x=5..7 rigid; 66x x=6..7 rigid.
- n=15: 77x x=7..8 rigid, pairs [0,5],[1,4],[2,3] (a = 7, sizes 5648/18416/8693), predicted rigid in the script, no full translator from collapse 67.
Did not re-run the 4-minute n=14 55x job; the committed output files agree with the claim text. The claim is reproduced, and extends to n=12,13 and to a=7 at n=15.

## True?
No counterexample found. Gaps:
- "Criterion" is a prediction rule checked on 4 values of a (+77x by me); a=4 is the only merged case, so the rule fits one positive datum. The "prediction made before running" is not checkable from the record (script contains the rule, not a timestamp).
- Formula k = 2x+3-a is stated for x range 5..8 (a=5), 6..8 (a=6); the statement "all orbits closed" holds in my runs too (no CAP).
- The claim that `aax@o -> aa(x-1)@(o+1)` holds for 44x, 55x, 66x is "see the table", but the table lists the drift only for the seed (aa(a+1)@(o-1)); the drift for general x is not shown separately. Indirectly supported by the pair structure at x>a.
- Mechanism for 34 (3333<->3403) is admitted unexplained; fine, but then the title's "closes exactly when" is a computed-criterion statement, not a theorem. The title says "exactly" and the body says "not proved": tone the title.
- Untested: a=2 (22x), a=8,9 (only 77 done, by me), n>=16 for any, and 34x/45x (correctly flagged silent).

## New?
grep of research/ (FINDINGS, HYPOTHESES, EXPERIMENTS) for collapse, translator, 55x, 66x, 44x: "collapse" hits are the unrelated end-pair collapse (F-022/F-029 region); 44x appears only as "not run" in a Limits line at research/EXPERIMENTS.md:21 and in E-065. Nothing on 55x/66x closure or the collapse-to-34 reason. Author's novelty claim stands. Literature hits are incidental.

## Evidenced?
Mostly. Specific n, x ranges, orbit sizes, offsets and the reproduction commands are given. Missing: n other than 14 for 55x/66x (I supply 12,13; 77x at 15), a record of the pre-registered prediction, and the null for "only a=4 special" (author defers to skeptic). The 10/24 negative result on the broader criterion is a useful stated control.

## Required for acceptance
1. Soften the title ("exactly") to match "computed; not proved", or add a=2 and a=8 cases.
2. Add the n=12, 13 (55x, 66x) and n=15 (77x) runs to the evidence; they agree.
3. State the drift for x>a in 44x/55x/66x explicitly (table or script output), not only the seed's neighbours.
4. Say how "prediction before running" is verifiable (commit order or the script printing PREDICT, which it does).
