# Experimentalist notebook (rewritten each round)

## What I now believe (after round 049)
- T10 (ii) control done: n = 7 c1, c2, 20 000 expansions each, reproduces E-149 (16 and 9 failing key-keepers, parent depth 7-8). Contracting the seeds (undirected shadow; walk length 15-18 follows from parent depth 7-8, not independent) and testing the explored graph: c2 all 9 pendant but all unexpanded (cap artifact, not dead ends; round-049 response); c1 14 of 16 in the seeds' biconnected block but only on walks of length 15-18, 2 pendant. Shortest merge walk in both graphs is 4. So no failing key-keeper is on a shortest merge path; the tally can see failures, hence the n = 8 c2 null (E-153) is a real null for depth <= 8.
- Earlier: 10 F-041 n=8 merges need no J != 0 step (E-148, E-154); guarded depth-4 BFS has 0 J != 0 edges; J != 0 key-preserving steps (E-145) at distance >= 5.
- Failing key-keepers are deep (parent depth 7-8); 13 of 25 children unexpanded; whether they leave the derived class is still open (E-152; toolsmith T10 i).

## What I tried
- 049: experimentalist_deepreplay.py (047 script + edge list in checkpoint), experimentalist_blocks.py (contract seeds, biconnected block, walk length with e removed). 2 slices x 2 classes in parallel, about 1 400 s wall.
- 047 deepreplay n=8 c2 depth 8; 045 t10b/_sweep; 043 tally; 041 keyoff; 038 d3table; 035 dhist; 033 bothdie; 030 BFS; 027 W walks.
- Lessons: background + `timeout N bash -c 'until grep ...'` poll; output to files, python -u; never `pkill -f` (it killed my own shell); a --budget-sec is needed under `timeout 10m` or the run dies with no checkpoint (lost 10 min on c1); --max-exp is in the script.

## What I would do next
1. Depth 9 for n = 7 c1, c2 (depth 9 is 54/85 percent undone for c1/c2; checkpoints in /tmp are gone next session; commit compressed edge lists if wanted; rerun ~25 min each, two slices); see whether new failing edges appear and whether they join the block.
2. Depth 9 n = 8 c2 overnight (proposal in the 047 submission).
3. Inspect the second E-094 rejection (path (17,8,5,6,8,8,2,5) vertex 5): still never examined.
4. Check my BFS against search.meetingPoints; n = 9 c0 shape tally.
- Watch: caps are not verdicts; the explored graph gives upper bounds on distances (a shorter route through a failing edge could exist in unexplored parts); taint is first-reach only.
