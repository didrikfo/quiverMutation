# Review of workshop/rounds/047/experimentalist.md

referee: skeptic · round: 047
verdict: minor revision

## Reproduction

Re-ran slice 1 only (the author's invited check), fresh checkpoint in the scratchpad: `timeout 590 .venv/bin/python -u workshop/rounds/047/experimentalist_deepreplay.py 8 --class 2 --depth 8 --budget-sec 450 --ckpt <new>`; 460 s wall (7m40s), exit 2. Depth 1-7 lines identical to the author's (expanded 18/96/290/734/1764/4108/9736, next 14580 at depth 7). The budget stop is wall-clock, so mine stopped at depth 8 pos 4938 (14 674 expansions) vs the author's 5089 (14 825). At that point: 1 rejection (path (6,1,1,4,3,1,7,4), vertex 4, the same as the author's first), merge 24 044, tree 36 984, guard-tilt 61 028, noguard-NOTtilt 1, FAILING list empty, GATE+TILT BUT KEY MOVES 0. Consistent with the author's slice-1 numbers (24 306 merge at 14 825), and with the stated one-failure state. Slice 2 (about 320 s) was not re-run, so the final 41 424 / 63 205 / 2 rest on the author's output file only; the partial-slice tallies are consistent with them.

## True?

No counterexample found. Two things the author did not state:

1. Tree-edge count does not reconcile with the distinct-algebra count. Tree edges + starts should equal distinct algebras. Full run: 63 205 + 18 = 63 223, reported 63 221 (off by 2). My slice: 36 984 + 18 = 37 002, reported 37 001 (off by 1). So the tree/merge split is off by a small amount, at least. Possible causes (start algebras that coincide, or a child counted "new" twice). The tilting tally (the claim) is unaffected, but the merge-edge figure 41 424 is not exact if the tree one is not. Compare E-094's unexplained +7.
2. The taint bookkeeping is vacuous here: with 0 failing guarded edges there is nothing to taint. It adds no evidence beyond the edge table; "tainted nodes: 0" should not be offered as a separate finding.

Not a flaw, but worth saying: the null is for key-kept edges, and the key guard is exactly what refuses the 2 failures, so "no merge path uses a non-tilting step" follows directly from "every kept edge is tilting". Fine.

## New?

Mostly E-094 (same walk, same 24 316 / 63 221 / 2). E-149 (16 of 80 978 n = 7 failing key-keepers at parent depth 7-8) and E-152 are the relevant contrast. New: only the merge/tree split and the explicit 0 failing key-kept edges over the full depth-8 expansion. Note E-094's 'guard-tilt 89 189' was over the first 20 899 expansions; the author's 104 629 is the whole walk; no conflict. Modest novelty, correctly labelled by the author.

## Evidenced?

Mostly. The edge table, the scope and the refutation condition are specific. Missing: the output file's slice order (slice 2 first) is stated; the final-state tally was produced across a resume and only count-agreement with E-094 checks it (author says so). The tree/merge reconciliation above is absent.

## Scope

Title and claim match what was checked (n = 8, class 2, depth 8, not closed). The "Consistency" paragraph says E-149's n = 8 showed none at 1 200 expansions; fine. The title's "no merge path uses a non-tilting step" is stronger than the tally: it is about the walk's edges, with first-reach taint, and the author caveats this. Narrowed wording: "no key-kept edge among the 104 629 generated in the depth-8 walk fails `tiltingPlus`".

## Required for acceptance

1. Reconcile tree edges + starts against 63 221 distinct (off by 2; mine off by 1 at 14 674 expansions), or state it as an unexplained discrepancy and say which of the split numbers it affects.
2. Drop or demote "Tainted nodes: 0" as evidence (vacuous when no failing guarded edges exist).
3. State that the final tally spans a resume and that slice 1 (re-run by the referee) matches the partial tally; the slice 2 portion was not independently re-run.
