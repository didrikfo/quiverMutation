# Theorist notebook (rewritten round 111)

## What I believe now
- T5 loop case. M1 (pivot-generator type, d f = ab) is proved by change of generators + Kunneth (E-290). M2 (free degree-1 loop, d f = 0) is proved in 111: filter by g = #{bbar, fbar, theta1}; d never raises g; gr has d0 theta1 = -v fbar, d0 theta2 = v u, a disjoint-alphabet Lemma A (E-240) with X1 = {v}, Y = {fbar, u}; finite filtrations per cell, so the q.i. lifts. Controls: planted wrong g caught, star ideal H != 0 (37/66 cells).
- M2 is a formal dg algebra, not a Construction-2.2 model of a finite-dim algebra (d f = 0 on a degree-1 letter). So M1/M2 together still do not decide H-023 (every gate-admitted J = 0 step is a derived equivalence); library steps never meet loops anyway (E-211, E-279).
- Degree-0 loop g (w = 3): the engine runs (E-284's 'not testable' only reflects the weight assertion), same dim I as M2, H(I) = 0 to weight 12, and the Lemma-A proof is degree-free, so the cancellation does not need Oppermann's 'no degree-0 loop'. What breaks is elsewhere (u has degree -1: kM not non-negatively graded); location unknown.
- Older: Hom control for M1 (T = P_j + cone, Cartan [[1,1],[1,2]], planted errors say NO); 106 six models cancel to weight 12; 103 loop of degree >= 1 at i forces an oriented cycle (E-279); Anick numbers UNVERIFIED (secondary).

## What I tried
- 111: letters/d of M2; g-filtration proof + machine check of d0 and true d (identical results); deg-0 loop parent; planted-filtration control.
- 109: M1 proof and Hom control; 106: loop filtration script (the engine all of this runs on: workshop/rounds/106/theorist_loopfilt.py).

## Next
- Same g-filtration for M4 (shared pivot) and M5 (product pivot ebar bbar); then classify loop models (pivot generator / Lemma-A leading part) and state one theorem.
- Find where non-negative grading enters after Prop. 7.6 (scholar), to say what a degree-0 loop breaks.
- T4: shape of the 14 width-5 middle algebras; 6 step-2 failures (w6/w7); one-step rule for width-5 class-B rules (E-266).
- Blind spots: dg formulas are the printed ones (E-202), not re-derived; Lemma A used for R = k^2 with a field Kunneth-free skeleton proof only (referee-accepted E-240); no Hom/End check for M2 (infinite-dimensional); deg-0 result is a single parent, one weight choice.
