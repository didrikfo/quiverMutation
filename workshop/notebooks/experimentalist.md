# Experimentalist notebook (rewritten each round)

## What I now believe (after round 043)
- E-143's shape (|out v| = 1, u = e_w - e_i, C'e_v = e_v) is an n = 6, 7 phenomenon: at n = 8 c0, c1 all 27 J != 0 steps found (520 s cap, depth 7-9) have |out v| = 2 (or |supp J| = 2); c2 had none. P1/P2 vacuous there.
- Q lowest term x^2 coeff 1: n6 146/146, n7 59/59, n8 c0 20/20; FAILS n8 c1 (x^3 in 5 of 7). Still never Q = 0 (key moves; E-140).
- n6: F order 8, s=1 -> c_2 = 0 (56/56), s=-2 -> c_2 = 1. n7 c0: F order 5, c_2 = 0 49/49. n7 c1: order 12, 8 of 10 steps have NO s with |s|<=8: the orbit relation is not universal.
- Earlier (041): key guard off, 0 of ~280 J != 0 steps keep the class key; random non-LNA parents do (55/1264). 038: n8 c0 J != 0 rows had out(i) = 3 (E-138, 192 steps).

## What I tried
- 043: experimentalist_tally.py (dump + reduce in one pass, --plan). Six jobs on 4 cores, 520 s each. Counts load-dependent.
- 041: keyoff, control. 038 d3table. 035 dhist. 033 bothdie. 030 BFS. 027 W walks.
- Lessons: output to files, `python -u`, never `pkill -f`; Bash waits >2 min need an until-loop with timeout 115, repeated.

## What I would do next
1. Print and save the 5 n8 c1 x^3 steps (relations, J, u); `off` mode untested: n8 c0/c1 guard off depth 6 as overnight (more steps than 27).
2. n = 9 c0 shape tally (|out v|, u).
3. Vet the 55 random control hits (relations, derived equivalence).
4. n7 c1 steps lacking s: what is e_i in F-orbit terms?
- Watch: caps are not verdicts; empty cells are vacuous; E-138 had 192 steps vs my 20 at n8 c0.
