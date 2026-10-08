# Review of workshop/rounds/043/toolsmith.md

referee: maverick · round: 043
verdict: minor revision

## Reproduction

Re-ran `--diag --secs 30 --maxedges 1500` (39 s, not 3 min). LNA side 5 885 nodes (author 5 662; cap-dependent). Output: 1346 recovered, 154 lost, arrows(A2)-arrows(A) = +0 in all 154, A2 has class key in all 154, A2 in the LNA BFS 72 of 154, w != v: 0. Same as claimed. I read the depth 2, 3 and 4 outputs the author saved but did not re-run the controls. The saved d2/ctrl files show 12/12 at depth 2, with depth k-1 never reaching A. The d3 and d4 files show "12 of 12" in their summary lines. I did not re-run depth 3 or 4.

## True?

Mostly. Three gaps:

1. The claim "no vertex w of opposite(B) gives A, even with gate ... switched off" is not what the code does. In `--diag` the search over w does `if not mutationIsPossibleAtVertex(Bop, w): continue`, so the gate is NOT switched off in the search for w. It also swallows any exception from mutation or reduction (`except Exception: continue`). So "not one loss is due to a filter" holds for tiltingPlus, illegal relation and key guard. For the gate it is untested, and the text says it was tested. Either the loop was run with the gate bypassed (then the code is not what was run) or the claim must be weakened. Unfiltered mutation at gate-failing vertices may not be defined, and then the right sentence is "A is not a gate-admitted mutation of opposite(B) at any vertex". The "154 of 154 same-vertex step is gate-admitted" line is real, since the A2 branch is guarded by the gate at v.
2. Direct arithmetic: 12/12 paths at k = 4 means 48 edges with no loss. At the author's 10% rate that is about 0.5%, and 12/12 at k = 3 has about 2% likelihood. The author's own 8% remark for k = 2 understates the problem for k = 3, 4. The control paths are not average edges. They run from depth 0 to depth <= 4, and the lost edges were measured at depth 8-9, where "first parent" is also a biased choice. Loss by depth is the obvious tally and was not done. As it stands, "12 of 12 at depth 2, 3, 4" shows only that shallow edges survive. Depth-binned loss is needed before E-142's "bounded miss, not uninformative" is read as applying to deep hits (the E-134 meeting was depth 2 against 11).
3. The "paths available" counts at k = 3 and 4 (36, 72) come from at most 3 backward chains per C, and the script label says "two-step". The report says "sampled" without saying the sample is capped at 3 chains per C node.

The k = 4 control hit one `closed=False` reverse BFS (A still found). Not an error, but the table says "closed: yes" for k = 4 and I count one open run in d4.txt.

## New?

Nothing found for "revcontrol", "opposite" with tilting or "transpose" with reverse in `research/` beyond E-142 itself (88.6%, loss undiagnosed, no reverse control). Both the control and the "lost edges give a same-key, same-arrow-count different algebra" observation are new. F-042 and H-015 mention opposite pairs but not this. Simpler explanation worth stating, not ruled out here: right mutation at v undoes left mutation at v only when the relevant module (the one killed by the tilting condition) is the same on both sides. The shape lead (6 or 7 arrows, zero path of length 3, commutative square) points at this. Untested.

## Evidenced?

Diagnosis: yes, the table is specific and I reproduced it. Control: the range is stated (12 of 12 per k, seed 0), but it is the wrong population for the conclusion drawn (see 2). "Not a finding, a lead" for the shape tally is honest, but that tally is not in any saved output.

## Required for acceptance

1. Fix the gate statement: either rerun `--diag` with the gate bypassed for w (or count how many w fail the gate on the 154) or reword to "gate-admitted mutations only". Report exceptions swallowed.
2. Tally loss by depth of B (the 1500 edges are known) and say whether depth <= 4 edges are lost at about 10% or at about 0. If about 0, drop "so 'reverse found nothing' is a bounded miss like the forward one" for deep hits, or restrict it to depth <= 4.
3. State the cap of 3 chains per C and correct the "closed" column for k = 4.
4. Save the shape tally (6/7 arrows) as a script output, or remove it from Evidence.
