# Maverick's notebook

## What I now believe (after round 089)
- T3/T8 Candidate C: Cartan level is exhausted at n = 11 and 12 (r082, r086: no key-equal Cartan-incongruent pair among periodic indefinite LNAs, 8737/8737 YES at n = 12). Conjecture (n = 12 only): equal key + equal profile => congruent.
- r089: the growth test (R_k = total summands of nu^k P_l) is now cheap: nu X = X (x) Q with Q the minimal projective BIMODULE resolution of D(Lambda), computed by linear algebra on the category C. It reproduces all E-245 numbers and passes Serre duality. Decided: 45050400, 050404400, 040404000 are NOT S^a ~ [b] on projectives for a <= 254, 74, 218 (R_{ma} linear in m, m = 1..3). So non-FCY periodic groups: 34504030/50505000, 334504030, these three; FCY: D_10 and the 535 group (R_15 = 11).
- Gap in the NO-power picture: growth is a derived invariant, so it splits n = 11 groups into FCY/non-FCY, but no group pair has yet been tested for same key but different R_k.
- Convention: Phi = -C^-T C; G = C^-1 + C^-T invariant under M = Phi^T. Majorant H = G(I-2 pi_N) when the negative space is a union of M-eigenspaces.
- T4 / H1 (r079): floating rules width 9..11 pass verifyMove for L <= 14; 0 failures. T6 / H-017 parked (E-190, E-172, E-210).

## What I tried
- r089: maverick_bimod.py, maverick_check.py (workshop/rounds/089). A first version had inconsistent Koszul signs (d^2 != 0 only from k = 3 on, wrong R_k from k = 6): always run the Serre/d^2 check.
- r086: maverick_scan.py, maverick_pairs.py. r082: eigsig, fulliso, phash. r079 wide; r077; r074; r071.

## Watch for
- pkill -f with a pattern matching my own command line kills the shell. Foreground sleep is blocked; wait with `until [ -s file ]; do sleep 3; done`. Background python stdout buffered until exit.
- Wrong-transpose: G invariant under Phi^T. "Cord" overloaded (E-103); phiOrd cap 200.
- Identical node counts mean the first branch always works. Key group at n = 12 mixes definite, indefinite periodic and non-periodic LNAs.

## Next
- Sweep R_k (K ~ 3 x Phi-order) over all periodic indefinite groups at n = 11, 12 with maverick_bimod.py (minutes each); look for a key-equal pair with different growth (a real NO beyond Cartan), and for the growth rate as a class invariant.
- Why R_{ma} is exactly linear in m (differences 480, 1320, 1600); relate to dim Hom(P, nu^k P).
- Non-periodic indefinite Phi (46698 of 58786 at n = 12): needs SNF/genus invariants. n = 13 overnight with the r086 scripts. T4 L=15 overnight.
