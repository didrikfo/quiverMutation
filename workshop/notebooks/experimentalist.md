# Experimentalist notebook (rewritten each round)

## What I now believe (after round 054)
- T10 (ii) group A is closed at replay level: 5 total-7 paths from [0,5,0,4,0,3,3,0] (05040330) to members of orbit 33460000 (e.g. 4 3 2 1 7 6 9, left steps negative), 35/35 edges gate, J = 0 (`tiltingPlus`), key-keeping, ends on the target row. All 5 meetings are non-labelled (phi nontrivial), so labelled `meetingPoints` misses them.
- Cheap: forward depth 3/4/5 = 131/387/1091 quivers (4/23/108 s); backward depth 3 for 42 targets 117 s. The guard prunes hard; the earlier depth-5 merges.py cost (80-150 s per member) was the all-122 sweep.
- The checkpoint "orbit" label 33460000 covers many relation rows (45055000, 55504400, 60504030, ...), 42 distinct members.
- Earlier (051): F-037 5 paths replayed 35/35; n = 10 merges.py depths 3-5 find no link (E-162). 049: n = 8 F-041 merges need no J != 0 step; n = 8 c2 depth 8 null real.

## What I tried
- 054: `experimentalist_amerge.py` (fwd/back reach stages, WL hash + VF2 iso join, signed-path inverse and replay) and `_join.py`.
- 051: merges.py 10 chunks, f037replay, a_meet. Earlier: deepreplay, blocks, t10b, tally, keyoff, d3table, dhist, bothdie, BFS, W walks.
- Lessons: background + poll with until-loops (no sleep chains); python -u; never pkill -f; an LNA end is recognised by lm.asRelLengths; a join up to relabelling needs invariant hash plus exact iso, and a transported inverse path (reverse, negate sign, rename by phi).

## What I would do next
1. Same join for F-037 (34504030 -> 50505000, 7 steps) with this script (start/target parametrised); find whether a 4 + 3 join exists, to cross-check.
2. Check other depth splits (3+4, dual starts) for more group-A paths, and whether any total-7 join uses a J != 0 step (it would show as J0 < edges); extend to total 8 using fwd 5.
3. Depth 9 for n = 7 c1, c2 (049); second E-094 rejection; n = 9 c0 shape tally.
- Watch: only the paths found, not all paths, were tested for J = 0; 7 shortest relies on E-033.
