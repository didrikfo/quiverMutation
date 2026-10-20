# Maverick's notebook

## What I now believe (after round 077)
- T6 / H-017: r077 ran 12 more defect-1 members (n = 9, L = 5) of the six left-over LNAs: 12/12 hit, 12/12 none at depth 4; 17/17 over 7 of 9 LNAs with E-203 (19/19 with referee's two; 4444400 has no cords >= 1 member). Controls test walk invertibility only; defect <= 0 detection untested. I recommended parking T6 at round 080 (no breadth slots, a scheduling decision; label: no counterexample at depth 5); reopen only for a proof route, the overnight depth-7 census finding a new cell, an independent defect <= 0 member, or a cheap n = 10 control.

(Earlier, after round 074:)
- T3/T8 Candidate C at n = 11 stays undecided. Cartan-level probes are exhausted in the cheap range (E-190, E-197, E-209): SNF(xC+C^T), Z-conjugacy of Phi, per-factor (rank, sig, det, theta). r074 adds explicit equivariant isometries: the definite pieces (ranks 3, 2) of all 15 members of the 15-group are isometric to member 0 (15/15 YES); rank-6 piece (indefinite (4,2), det 7) and the glue NOT run.
- Convention: Phi = -C^-T C; G = C^-1 + C^-T invariant under M = Phi^T; Cartan congruence = M-equivariant isometry of G. A NO there proves derived inequivalence; a YES proves nothing.
- Stage 1 algorithm (maverick_isometry.py): generators v_i of min norm spanning Q^n under M, enumerate images by norm + Gram constraints, P = W(u) W(v)^-1. Exact, fast (n=10 in <1 s). n = 10 definite periodic LNAs: only two key groups (256 A_10-type, 321 D_10-type), each a single isometry class (575/575 YES); 90000000 ~ 00000030 gets an explicit P.
- Test has no demonstrated NO-power beyond signature: no known inequivalent n=10/11 pair of equal key+signature. Artificial NO ((G,M) vs (G,-M)) fires; Galois twists M^k are YES.
- Real size of C: orbit+mirror classes 15-group 9, 24-group 10, 64-group 25, 535-group 213. The 16 unplaced n=11 LNAs lie in none of the four groups.
- Still true: n=10 Phi^18 group has 2 classes (E-190); HH^* = k (E-172); T6 controls complete, parked; E-210: 90000000 / 00000030 deletions all land in one n = 9 class.

## What I tried
- r077: maverick_defect1_rest.py (shard of 069 script, env LO); 4 parallel jobs of 2 members each fit in 10 min (4 cores, 250-530 s per member).
- r074: maverick_isometry, _survey10, _negctl, _pieces15 (all < 1 min). r071: snfx, orbits, conj, pieces, place16, control10. Earlier: defect-1 control, finer/factors/qform, cartanC, n11periodic, cordcontrol, group, ctrlpairs.

## Watch for
- pkill -f with a pattern matching my own command line kills the shell: do not use it.
- Wrong-transpose: G is invariant under Phi^T (assert in data()).
- "non-Cartan" test names: isometry stages are Cartan-level in disguise; only a J = 0 join search is a different kind of evidence.
- grep research/ before announcing; "cord" overloaded (E-103); phiOrd cap 200.

## Next
- Stage 2: rank-6 piece as rank-1 hermitian lattice over Z[zeta_7] (h = 1; unit-norm criterion in the totally-positive units of Q(zeta_7)^+ mod squares); ask theorist/scholar to confirm first. Then stage 3 glue (discriminant groups).
- Toolsmith's J = 0 join search with 90000000 ~ 00000030 as control remains the only non-Cartan route.
- Overnight: forward census depth 7 (T6); S-1 n = 15 K0 = 5.
