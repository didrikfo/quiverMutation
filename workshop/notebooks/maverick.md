# Maverick's notebook

## What I now believe (after round 082)
- T3/T8 Candidate C at n = 11: Cartan level is exhausted. r082: 8 of 8 non-identity orbit+mirror reps (plus identity) of the 15-group are Cartan-congruent to rep 0 on the FULL indefinite lattice (explicit P, exact check, P C_a P^T = C_b); per-eigenvalue signature is one profile for all 15. Only a non-Cartan separator or a join can move it. The tests have no demonstrated NO-power beyond signature (only artificial (G,M^5) NOs).
- Convention: Phi = -C^-T C; G = C^-1 + C^-T invariant under M = Phi^T. Majorant H = G(I-2 pi_N) works when the negative space is a union of M-eigenspaces (true here, angle 1/14).
- T4 / H1 (r079): floating rules width 9..11 pass verifyMove for L <= 14; L=15 partial; 0 failures. Same engine as E-159..E-207, no independence.
- T6 / H-017: parked (r077): 17/17 defect-1 members hit; reopen only for a proof route, a new census cell, or an n=10 control.
- Still true: n=10 Phi^18 group has 2 classes (E-190); HH^* = k (E-172); E-210 deletions land in one n = 9 class.

## What I tried
- r082 response: maverick_phash.py -> maverick_P.txt (node cap 2e6 not hit, Bmax 2.4583, sha of all nine P).
- r082: maverick_eigsig.py, maverick_fulliso.py, maverick_fulliso_ctl.py (workshop/rounds/082). 
- r079 maverick_wide.py; r077 defect1_rest; r074 isometry/survey10/negctl/pieces15; r071 snfx, orbits, conj, pieces, place16, control10.

## Watch for
- pkill -f with a pattern matching my own command line kills the shell. `sleep` chains are blocked.
- Wrong-transpose: G invariant under Phi^T. "Cord" overloaded (E-103); phiOrd cap 200.
- M^k twists must use k coprime to the full order (42 at n=11, 18 at n=10), else angles collide.
- Identical node counts across targets (138) mean the first branch always works: not evidence of power.

## Next
- Look for real power: a pair of LNAs with equal key, equal eigenvalue-signature profile, and Cartan-incongruent (scan n = 12, 13 key groups with fulliso). If none exists, say so as a theorem-shaped conjecture.
- Toolsmith J = 0 join search with 90000000 ~ 00000030 as control (still the only non-Cartan stage).
- T4 L=15 widths 9..11 overnight; ask whether the 18 rules per width are a parametric family.
