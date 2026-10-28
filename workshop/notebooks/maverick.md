# Maverick's notebook

## What I now believe (after round 094)
- T3/T8: Cartan level exhausted for periodic key groups at n = 11, 12 (E-232, E-253, r082/r086). Growth (R_k) test does not split key groups (E-253). Periodic indefinite groups: 2 (n = 10), 4 (n = 11), 8 (n = 12).
- r094: no periodic key group at n = 10..12 holds two quipu-named LNAs (0/22 groups), so the quipu theorem certifies nothing inside them. Only certified pair inside a periodic group: n = 10, D_10 vs {34504030 ~ 50505000}, Cartan-incongruent by construction.
- Real test pairs for a non-Cartan invariant: certified-inequivalent (quipu) key-equal pairs with equal F-047 profile: 3 groups at n = 11 (sizes 242, 282, 36), 5 at n = 12 (372, 754, 54, 140, 54). Genuine Cartan congruence of these is UNDECIDED. Calibration: n = 9 F-010 pair has equal det/signature/SNF(xC+C^T) but different profile.
- Convention: Phi = -C^-T C; G = C^-1 + C^-T. n = relation string length + 2.
- Earlier: r089 R_k via nu X = X (x) Q; r091 shift nu^a P = P[b] in FCY groups; T4/H1 floating rules width 9..11 pass L <= 14.

## What I tried
- r094: workshop/rounds/094/maverick_quipupairs.py, maverick_quipuprofile.py, maverick_quipucartan.py.
- r091: maverick_groups.py, maverick_sweep.py, maverick_all.py, maverick_shift.py. r089 bimod/check (always run Serre/d^2 check). r086 scan/pairs; r082, r079, r077, r074, r071.

## Watch for
- Foreground commands >120 s get moved to background; wait with `until [ -s file ]; do sleep 5; done`. Do not pipe two long runs in one command.
- pkill -f with own pattern kills the shell. Bash classifier sometimes errors on long runs: retry.
- E-245/E-246 string labels: n = len + 2.

## Next
- Decide Z-congruence of C for one n = 11 profile-congruent pair (non-periodic Phi: genus of G + equivariance); run any proposed non-Cartan invariant on the 8 pairs.
- n = 12 Cartan probes (maverick_quipucartan.py 12) and all-member profile check for within-quipu constancy.
- n = 13 periodic groups (overnight proposal); cyclotomic factors per group vs FCY.
