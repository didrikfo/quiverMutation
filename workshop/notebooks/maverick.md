# Maverick's notebook

## What I now believe (after round 086)
- T3/T8 Candidate C: Cartan level is exhausted at n = 11 (r082, 15-group, 8/8 YES) AND at n = 12 (r086): 8 key groups with indefinite periodic members (8745 LNAs), one signature profile each, 8737/8737 YES against the group rep, node cap never hit. No key-equal Cartan-incongruent pair exists among periodic indefinite LNAs at n <= 12, so the Cartan test has no demonstrated NO power beyond signature/key. Conjecture (n = 12 only): equal key + equal profile => congruent.
- n = 13 not run: key table ~6 min + classification + iso stage > 10 min; proposed overnight.
- Convention: Phi = -C^-T C; G = C^-1 + C^-T invariant under M = Phi^T. Majorant H = G(I-2 pi_N) works when the negative space is a union of M-eigenspaces (finite-order Phi).
- T4 / H1 (r079): floating rules width 9..11 pass verifyMove for L <= 14; L=15 partial; 0 failures. Same engine as before, no independence.
- T6 / H-017: parked (r077). Still true: n=10 Phi^18 group has 2 classes (E-190); HH^* = k (E-172); E-210 deletions land in one n = 9 class.

## What I tried
- r086: maverick_scan.py, maverick_pairs.py (workshop/rounds/086; outputs maverick_pairs12_all.txt).
- r082: maverick_eigsig.py, maverick_fulliso.py, maverick_fulliso_ctl.py, maverick_phash.py.
- r079 wide; r077 defect1_rest; r074 isometry/survey10/negctl/pieces15; r071 snfx, orbits, conj, pieces, place16, control10.

## Watch for
- pkill -f with a pattern matching my own command line kills the shell. Foreground `sleep` is blocked; wait with `until [ -s file ]; do sleep 3; done`. Background python stdout is buffered until exit.
- Wrong-transpose: G invariant under Phi^T. "Cord" overloaded (E-103); phiOrd cap 200.
- M^k twists need k coprime to the full order. Identical node counts mean the first branch always works.
- The key group at n = 12 mixes definite, indefinite periodic and non-periodic LNAs; only 14.9% of LNAs are indefinite periodic.

## Next
- Non-periodic indefinite Phi (46698 of 58786 at n = 12): is there a Cartan-incongruent key-equal pair there? Needs a majorant without finite order (or SNF/genus invariants).
- Toolsmith J = 0 join search with 90000000 ~ 00000030 as control (the only non-Cartan stage).
- n = 13 overnight with the two r086 scripts. T4 L=15 widths 9..11 overnight.
