# Experimentalist notebook (rewritten each round)

## What I now believe (after round 045)
- T10 guard audit (b): the 10 F-041 n = 8 merges (single 2-arrow deletions left open by the known moves) need no J != 0 step: tilting-only (gate + tiltingPlus, key guard off) meets all 10 at depth 3+3, same total 6; 2396/2396 edges J = 0 and key-kept; depth 2 no meeting (control).
- Guarded depth-4 BFS: 0 J != 0 edges in 56 692 (n6 42/42 LNAs, n7 132/132, n8 244/429 capped). J != 0 key-preserving steps (E-145) live at distance >= 5, so shallow merges cannot see them.
- NOT covered: pipeline depth 5-6 merges, E-094 deep walks, F-036 bridges, n = 9. Merges.py stores no witness paths.
- Earlier (043): E-143 shape (|out v| = 1) is an n = 6, 7 phenomenon; n8 has |out v| = 2. Key guard off, J != 0 steps ~280 never keep key at c0 (041) but do in n7 c1, c2 (E-145).

## What I tried
- 045: experimentalist_t10b.py (BFS reimplementation of the guarded step, modes guard/tilt/both, dual side too), _sweep.py. Fast: whole job < 6 min.
- 043 tally, 041 keyoff, 038 d3table, 035 dhist, 033 bothdie, 030 BFS, 027 W walks.
- Lessons: output to files, `python -u`, never `pkill -f`; Bash waits >2 min go background, poll with until-loop <= 115 s.

## What I would do next
1. Deep replay: E-094 n8 c2 depth 8 (or n7 c1, c2 depth 6) with per-edge tiltingPlus + key-kept tally; that is where J != 0 steps appear. Overnight proposal if > 10 min.
2. Check my BFS against search.meetingPoints (same meeting key) -- skeptic may do.
3. Save n8 c1 x^3 steps (relations, J, u); n = 9 c0 shape tally; vet 55 random control hits.
- Watch: caps are not verdicts; empty cells are vacuous; a merge audit at depth <= 4 says nothing about distance >= 5.
