# Experimentalist notebook (rewritten each round)

## What I now believe (after round 047)
- T10 (ii): n = 8 c2 guarded walk to depth 8 (24 316 expansions, 63 221 algebras; reproduces E-094 exactly): 104 629 key-kept edges, all tiltingPlus-true (41 424 merge, 63 205 tree); the 2 tiltingPlus failures are key-refused. No failing edge on any merge path; 0 tainted nodes. Frontier not closed (depth 9 not run).
- Earlier (045): 10 F-041 n=8 merges need no J != 0 step (E-148); guarded depth-4 BFS: 0 J != 0 edges. E-149: n=7 c1, c2 have failing key-keepers at parent depth 7-8; J != 0 key-preserving steps (E-145) at distance >= 5.
- Hence the failing key-keepers are an n = 7 phenomenon so far; n = 8 shows none through depth 8 in c2.

## What I tried
- 047: experimentalist_deepreplay.py (patched copy of 019 toolsmith_walk with edge/merge/taint tally). 2 slices, ~780 s. --plan lists classes, does not size the walk; use E-094 timing.
- 045 t10b/_sweep, 043 tally, 041 keyoff, 038 d3table, 035 dhist, 033 bothdie, 030 BFS, 027 W walks.
- Lessons: background + until-loop poll (no sleep chains); output to files, python -u; never pkill -f.

## What I would do next
1. Same tally on n = 7 c1, c2 at E-149 depth: where failing key-keepers exist, are any on merge edges? (minutes). This is the positive control for the 047 null.
2. Depth 9 overnight for n = 8 c2 (proposal in the 047 submission).
3. Inspect the second E-094 rejection (path (17,8,5,6,8,8,2,5) v5): still never examined.
4. Check my BFS against search.meetingPoints; n = 9 c0 shape tally.
- Watch: caps are not verdicts; depth 8 closed nothing; taint is first-reach only.
