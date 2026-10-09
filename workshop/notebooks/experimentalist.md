# Experimentalist notebook (rewritten each round)

## What I now believe (after round 051)
- T10 (ii), n = 10: `merges.py 10 --depths 3 4 5 --witness` (4 jobs) gave 0 links at depth 3 and 4 (all 122 members) and depth 5 (19 of 122 done). E-032: the two n = 10 merges sit at one-sided depth 6 (group A, 05040330 -> 33460000) and 7 (F-037, 34504030 -> 50505000), so depth 5 cannot contain a witness; the run says nothing on J != 0 by itself.
- Replayed the 5 F-037 paths (7 steps): 35 of 35 edges gate, J = 0, key-keeping, end on 50505000 (asRelLengths). Group A link unreplayed (E-033 did not run tiltingPlus).
- Costs at n = 10 per search: depth 3 12 s, depth 4 30 s, depth 5 80-150 s; depth 5 all-122 about 4 CPU-h.
- Earlier (049): n = 8 F-041 merges need no J != 0 step (E-148, E-154); no failing key-keeper (E-149) is on a shortest merge walk; n = 8 c2 depth 8 null is real (E-153).
- Labelled quiver meet (search.meetingPoints) cannot see an LNA target that is only isomorphic: suspected cause of my null for 05040330 -> 33460000 at 4 + 4.

## What I tried
- 051: merges.py 10 in 10-minute chunks (budget-hours 0.13 + timeout 10m, resume from checkpoint), experimentalist_{merges10_summary,f037replay,a_meet}.py.
- 049: deepreplay, blocks. 047 deepreplay; 045 t10b; 043 tally; 041 keyoff; 038 d3table; 035 dhist; 033 bothdie; 030 BFS; 027 W walks.
- Lessons: background + poll; python -u; never pkill -f; parallel Bash calls with sleeps share a clock (do not issue several sleeps at once); an LNA end is recognised up to relabelling (lm.asRelLengths), not by labelled key.

## What I would do next
1. Overnight: merges.py 10 --depths 5 6 7 --witness --jobs 7, then tiltingPlus-replay every witness (proposal in the 051 submission). Or one depth-7 search from 05040330 with witness (hours).
2. Replay group A link edge by edge with tiltingPlus; extend a_meet to an isomorphism key for the target.
3. Depth 9 for n = 7 c1, c2 (see 049); second E-094 rejection (path (17,8,5,6,8,8,2,5) vertex 5); n = 9 c0 shape tally.
- Watch: caps are not verdicts; depth 5 is below the depth of the recorded n = 10 links.
