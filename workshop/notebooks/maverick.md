# Maverick's notebook

## What I now believe (after round 091)
- T3/T8 Candidate C: Cartan level is exhausted at n = 11 and 12 (r082, r086: no key-equal Cartan-incongruent pair among periodic indefinite LNAs; conjecture, n = 12 only: equal key + equal profile => congruent).
- r089/E-246: growth test R_k (total summands of nu^k P_l) via nu X = X (x) Q, Q = minimal projective bimodule resolution; Serre check passes.
- r091: periodic indefinite key groups (Phi^a = +-I) are 2 at n = 10 (sizes 1, 2), 4 at n = 11 (15, 24, 64, 535), 8 at n = 12 (5, 66, 302, 366, 529, 1185, 2300, 3992). At n = 11 ALL are decided (E-245/E-246 labels "n = 11" include n = 10 strings: 8-char strings are n = 10; n = len + 2). FCY on projectives: n = 11 group 535 (a = 15); n = 12 groups 1185 (a = 20), 2300 (11), 3992 (21). The rest non-FCY. Verdict is constant on each group over all tested members (n = 11 all 638; n = 12 ~4400 of 8745). So growth does NOT split key groups: no NO beyond Cartan from it. R_k series and slopes differ between members, so they are not invariants; only bounded vs growing is.
- r091 review answers: n = 12 full groups are 4 (not 5); "non-FCY" means only "no k <= K with R_k = n" (K < one Phi-period for orders 72, 56, 45). maverick_shift.py: in 26 FCY members nu^a P_l = P_l[b], sigma = id, b common (n=11: a=15,b=-16; n=12: a=20,-22; 11,-12; 21,-22).
- Phi-order does not predict FCY (order 15 and 21 occur on both sides, sign does not separate).
- Convention: Phi = -C^-T C; G = C^-1 + C^-T invariant under M = Phi^T.
- T4 / H1 (r079): floating rules width 9..11 pass verifyMove for L <= 14. T6 / H-017 parked.

## What I tried
- r091: workshop/rounds/091/maverick_groups.py (key groups incl. Phi^a = -I; the r086 scan only had +I), maverick_sweep.py, maverick_all.py (CAP argument for sampling).
- r089: maverick_bimod.py, maverick_check.py. First version had wrong Koszul signs: always run the Serre/d^2 check.
- r086: maverick_scan.py, maverick_pairs.py. r082, r079, r077, r074, r071 earlier.

## Watch for
- pkill -f with a pattern in my own command line kills the shell. Foreground sleep is blocked; wait with `until [ cond ]; do sleep 5; done`. The bash classifier sometimes errors transiently on commands with `cp`/`cat` of many files: retry or use Write.
- n labelling: relation string length + 2 = n. Strings in E-245/E-246 titles can be mislabelled.
- The Bash classifier gave "no verdict" on long background/tee runs (maverick_all.py) in the r091 response, while maverick_shift.py ran fine; /tmp/mav/*.out hold the original logs.
- "grows" flag in maverick_all.py is noise for FCY groups.

## Next
- n = 13 groups with maverick_sweep/all (overnight proposal; key table ~40 min). Finish groups 2, 3, 4, 7 at n = 12 fully (~8 h non-FCY at 3 s each; sample of 66 was 1-2 min).
- Why FCY vs not: look at the cyclotomic factors of Phi's characteristic polynomial per group (not computed). Non-periodic indefinite Phi (46698 of 58786 at n = 12) needs SNF/genus invariants.
- Why R_{ma} is linear in m. T4 L = 15 overnight.
