# Maverick's notebook

## What I now believe (after round 071)
- T3/T8 Candidate C at n = 11 stays undecided, and Cartan-level probes are exhausted in the cheap range: SNF(xC+C^T) (x = -5..5; = SNF(xI - Phi^-1), a function of the Z[t]-module only), Z-conjugacy of Phi (14/14 explicit P in the 15-group), per-factor isometry class of (ker f(Phi^T), G) (rank, sig, det, theta): ONE value per group, all four groups (E-190 groups 535/64/24/15).
- Convention: Phi = -C^-T C; G = C^-1 + C^-T is invariant under Phi^T (not Phi); C^-1 = (I-Phi)^-1 G, so Cartan congruence = Phi^T-equivariant isometry of G.
- Real size of C: orbit+mirror classes (freeMoves.derivedOrbits free=True + mirrorRow): 15-group 9, 24-group 10, 64-group 25, 535-group 213. No orbit leaves its group.
- The 16 unplaced n = 11 LNAs (E-162): one quipu key, signature (9,1,1), Phi not periodic: in none of the four groups.
- Still true: n=10 Phi^18 group has 2 classes (E-190); HH^* = k (E-172); T6/H-017 controls complete (E-185, E-203), parked.
- r058: equivalent control joins at 4+4, ~3x per level; r050: 16 unplaced at n = 11.

## What I tried
- r071: maverick_snfx, _orbits, _conj, _pieces, _place16 (all < 70 s). r069 defect-1 control; r066 finer/factors/qform; r063 cartanC, n11periodic; earlier cordcontrol, group, ctrlpairs.

## Watch for
- Wrong-transpose mistakes: check numerically which of Phi, Phi^T preserves G before building a lattice invariant. My first piece run used Phi and gave 9 values on the 15-group (count = 9 orbit classes: unexamined coincidence).
- Power: no known inequivalent pair at n = 11 with equal key + signature + these invariants; a NO-capable control does not exist yet. n = 10 control 90000000 ~ 00000030 needs the indefinite filter dropped (D_10 G is definite).
- grep research/ and research/literature/ before announcing a dead end; "cord" is overloaded (E-103); phiOrd cap 200.

## Next
- Only non-Cartan separators can move C: J = 0 join between two of the 9 orbit+mirror classes of the 15-group (toolsmith; --plan first; n = 10 control first).
- Optional: staged explicit isometry on the 15-group (definite pieces finite, glue, indefinite rank-6 piece via Z[zeta_7]-units) for an E-154-style P between different orbit classes; or ask a theorist for the class-number statement.
- Examine whether the wrong-Phi tuple partition equals the orbit+mirror partition (cheap, may be a labelling invariant).
- Overnight: forward census depth 7 (T6); S-1 n = 15 K0 = 5.
