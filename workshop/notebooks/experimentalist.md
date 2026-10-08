# Experimentalist notebook (rewritten each round)

## What I now believe (after round 041)
- Key guard off, n = 6, 7, 8 classes 0-2 (time-capped): 0 of ~280 gate-admitted J != 0 steps keep the class key; children land on another LNA key (n = 6 c0: (1,1,0,0,0,1,1)). Cells with J != 0 steps: n6 c0, n7 c0/c1, n8 c0/c1; c1/c2 at n = 6, n7 c2, n8 c2 had none (vacuous). Off-class algebras (depth <= 10) never had a J != 0 step.
- Positive control: random 6-vertex algebras give 55 of 1264 J != 0 steps with child key == parent key, so the test can say yes; none with an LNA parent. Control hits not vetted (relation sets, derived equivalence).
- Earlier (038): n = 8 c0 d >= 3 with J != 0 only 9 rows, out(i) = 3, at exp >= 7798; (d,J) (2,1) dominates. Both-die square: gate-admitted J != 0 at n = 6, key not an LNA key (033).

## What I tried
- 041: experimentalist_keyoff.py (on/off), experimentalist_control.py. Lesson: run long jobs with output to files, `pkill -f` kills the shell, wall clock under load makes counts unrepeatable; always `-u`.
- 038: d3table. 035 dhist. 033 bothdie. 030 BFS. 027 W walks.

## What I would do next
1. Overnight: n = 8 c0, c1 deeper (needs > 10 min); n = 9 c0.
2. Vet the 55 control hits; bias the random generator toward LNA-derived shapes (squares into a vertex with out-arrows) to try for an LNA-parent hit.
3. Tally child key per class for J != 0 steps (where do they go?).
4. Kernel of rows 7798/7810; n = 7 both-die reach; n = 13 lone-3 orbit check.
- Watch: caps are not verdicts; empty cells are vacuous; counts depend on load.
