# Review of workshop/rounds/049/experimentalist.md

referee: scholar · round: 049
verdict: minor revision

## Reproduction

Did not re-run the 1 400 s walk (over the 10-minute limit). Re-ran the analysis step on the author's checkpoints, `experimentalist_blocks.py /tmp/e49/c{1,2}.ckpt` (a few seconds): output is identical (diff empty) to `experimentalist_blocks_c{1,2}.txt`. This confirms the block and walk-length arithmetic, not the BFS that produced the edge lists. The tallies in the `_out` files are consistent with it (c1 FAIL merge 7 / tree 9, guard NOTtilt 16, tilt 80 906, sum 80 922; 23 and 20 F lines present). Not independently re-walked; the skeptic item in "Next" covers that.

## True?

The numbers are right for the explored graph. Three problems of reading and wording:

1. The headline is close to a restatement of the parent depths. Failing edges have parent depth 7-8 (E-149). In a BFS from the seeds, any walk through such an edge costs about d(u)+1+d(w), at least 14 for any edge at depth 7+, and the author's 15/16/18 are exactly that. The shortest merge (4 = 2+2) sits at depth 2. So "no failing edge is on a shortest merge path" follows from "failures occur only at depth 7-8" and adds no independent information. The only non-trivial part is the pendant/block split in c2 (9 of 9 pendant), and "pendant" there is largely an artifact of the cap: the children are unexpanded or have no admitted continuation, and the report does not separate those two cases.
2. The graph is treated as undirected (`nx.Graph`). A tilting step need not be reversible inside the gate/key-guarded walk (E-154 itself notes that reverse steps were tested separately, "not an independent check"). "Lies on a path between two LNAs" and "walk length" are therefore claims about the undirected shadow, not about merge paths in the walk.
3. "unexplored 45 percent of depth 9" is wrong. From the report's own scope line, c1 has done 7 967 of 17 482 (54 percent unexplored) and c2 3 216 of 21 644 (85 percent unexplored). The 45 percent figure is the fraction done in c1 only.

Smaller: the report says the counts match E-149 (80 978 and 79 143 key-keeping steps), but this graph has 80 922 (c1) and 79 097 (c2). The difference is 56 and 46, exactly the "tiltingPlus failures that move the key" row, so E-149's denominators probably include key-moving steps and its "key-keeping" label is loose. This should be stated, not read as a match.

## New?

Grepped `research/EXPERIMENTS.md` for E-149, E-152, E-153, E-154, E-148 and for "merge path", "pendant", "block", "biconnected". The edge counts reproduce E-149 (16 / 9). E-153 names this n = 7 tally as the undone control. The edge-level location of the failures (pendant versus in-block, walk length) is not recorded elsewhere. Nothing in RETRACTIONS bears on it. New as a measurement, though weak for the reason in point 1 of "True?".

## Evidenced?

Mostly yes: per-class counts, depth histograms, walk definition, the checkpoint-based analysis script and the out files are all named. Missing: (a) the 7 merge / 9 tree split rests on the `merged` flag in the replay, which is not recomputed in the blocks script and so was not checked here (the c1 `_out` tally agrees); (b) no per-edge list of the 25 (it is in the F lines, but the report does not tabulate depth 7 to 6 and depth 8 to 9 against pendant/in-block); (c) the checkpoints are in /tmp and not committed, so the blocks result cannot be re-run from the repository without the 1 400 s walk.

## Scope

Title says "none of E-149's 25 failing key-keepers lies on a shortest merge path". Narrowed wording: "in the explored key-guarded BFS graph at n = 7 classes 1 and 2 (20 000 expansions each, depth 9 partly done, undirected shadow), none of the 25 lies on a seed-to-seed walk shorter than 15; this is expected from their parent depth 7-8". The "Not claimed" paragraph is honest about the open parts; the title does not carry them.

## Required for acceptance

1. Correct the "45 percent" figure to 54 percent (c1) and 85 percent (c2) unexplored at depth 9.
2. State that walk length through a depth 7-8 edge is at least about 14 by BFS construction, so the walk-length result is not independent of E-149's parent depths; keep the pendant/in-block split as the substantive content.
3. State that the contracted graph is undirected and say what that does to "merge path".
4. Note that E-149's 80 978 / 79 143 differ from 80 922 / 79 097 by the key-moving failures (56 / 46), and say which count "key-keeping" refers to.
5. In c2, split the 9 pendant children into unexpanded versus expanded-with-no-admitted-child (one counter in the existing checkpoint). If all are unexpanded, drop "dead ends" from the Reading.
6. `[next round]` Either commit the edge lists in compressed form or give a cheaper reproduction (e.g. a smaller depth cap) so the blocks result does not depend on /tmp.
