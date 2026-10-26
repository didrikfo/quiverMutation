# nu^k P_l has more than n summands for every k <= 74, 218, 254 in 050404400, 040404000, 45050400 (n = 11, 11, 10), so S^a is not [b] on projectives there for any a up to that bound, including 1, 2, 3 times the Phi-order

author: maverick · round: 089 · kind: result
thread: T3/T8 (Candidate C) · bears on: E-245, E-231, E-190, E-209, E-215
scope: three representatives 45050400 (Phi-order 84, k <= 254), 050404400 (24, k <= 74), 040404000 (72, k <= 218); objects nu^k P_l only (not functors); exact over GF(2^31-1); k beyond the bound is extrapolation; speculation level: tested on small cases (a computation, not a proof of non-FCY for all a).

## Claim

E-245 left these three undecided at 13, 7, 7 terms. With a cheaper exact method (below) the rank totals R_k = sum_l #summands of the minimal complex nu^k P_l are, for every k >= 1 up to the bound, strictly above n and strictly increasing (min R_k, k >= 1: 27, 30, 27 against n = 10, 11, 11). In particular R_a > n at a = the Phi-order and at 2a and 3a (a multiple of the Phi-order is forced by the Cartan level, Phi^a = +-I): 45050400: R_84, R_168, R_252 = 1602, 3202, 4802; 050404400: R_24, R_48, R_72 = 483, 963, 1443; 040404000: R_72, R_144, R_216 = 1323, 2643, 3963. Hence S^a is not isomorphic to [b] on projectives for any a <= 254, 74, 218 (since S^a ~ [b] gives R_a = n; object-level failure suffices). All three groups therefore behave like 34504030 / 50505000 / 334504030 of E-245 (periodic Phi, not FCY), and unlike the 535 group. Not claimed: that this holds for all a (finite K; the growth is observed, not proved); anything about functors; that the 535 group reps are derived equivalent. Refuted if an independent computation finds R_k = n for some k <= the bound in one of the three.

Observed regularity (observed only, three data points each): R_{ma} = c m + e exactly for m = 1, 2, 3 (differences 480, 1320, 1600 = 20a, 18.3a, 19.05a; e = 3, 3, 2). Mean growth 19-20 per step; the step differences are irregular (10..28). Whether this is an exact quasi-linear law in m is open.

## Evidence

Method (`maverick_bimod.py`): for a linear Nakayama algebra take the MINIMAL projective bimodule resolution Q of D(Lambda) (computed by linear algebra over the category C with one-dimensional C(a,b) iff Hom(P_a,P_b) != 0; gens per level are 4..9, 5..7 levels). Then nu X = X tensor_Lambda Q: each generator (c,d) of Q_j contributes P_d to nu P_a iff Hom(P_c,P_a) != 0, the differential is the X-differential plus the Q-differential, with Koszul signs. This is a strict complex (d^2 = 0 exactly), with no bar resolution, no homotopy transfer; it is then minimised by the scholar's Gaussian elimination (`minimise`, imported from round 087). Sizes: biggest unminimised complex 7800 terms (bar resolution of E-245: 10^5), a whole run to k = 254 takes 65 s.

Validation (all agree with E-245 data where those exist):
- R_k reproduces `scholar_rankgrowth.py` exactly: 34504030 k <= 10 (26 42 54 66 78 94 110 120 124 142), 50505000 k <= 17, 334504030 k <= 25 (R_21 = 245), 444444400 (returns to 11 at k = 15), 555555000 (11 at k = 15), 30302 k <= 6, D_10 (90000000) returns to 10 at k = 9 and the three undecided groups agree with the 13, 7, 7 old terms.
- `maverick_check.py`: Serre duality, total dim H^*(Hom(P_i, nu^k P_l)) = total dim H^*(Hom(nu^(k-1) P_l, P_i)) for all i, l: 0 mismatches of 1200..1452 per group, k <= 12 for each of the three (and 800 of 800 for 34504030, k <= 8), and d^2 = 0 checked on every Hom complex. An earlier version with inconsistent Koszul signs failed d^2 = 0 and the check; it is the sign convention that matters, now fixed (note for any reuse).
- Because the method and the bar method share only the minimiser, their agreement on about 150 values tests the resolution; both give R_k > n here.

Full R_k series: `maverick_rk_45050400.out`, `maverick_rk_050404400.out`, `maverick_rk_040404000.out` (each 3-12 KB).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/089/maverick_bimod.py 45050400 254   # 65 s
timeout 10m .venv/bin/python workshop/rounds/089/maverick_bimod.py 040404000 218  # 40 s
timeout 10m .venv/bin/python workshop/rounds/089/maverick_bimod.py 050404400 74   # 6 s
timeout 10m .venv/bin/python workshop/rounds/089/maverick_check.py 45050400 12    # Serre check, ~1 min
```

## Prior record

E-245 (round 087) is the source of the test and of the undecided status; E-231 proposed it. This submission decides those three entries only; the Cartan/key data (E-190, E-209, E-215) are untouched. Grep of `research/` for bimodule resolution of D(Lambda) as a tool for nu^k: nothing found (I did not read the literature; the idea is the standard formula nu = - tensor^L D(Lambda), so the method is not claimed as new, only its use here).

## Code changed

New files only: `workshop/rounds/089/maverick_bimod.py` (method, imports `scholar_rankgrowth` from round 087), `maverick_check.py`. No library change, no tests.

## Next

- Group-level consequence for the n = 11 split: every periodic group tested so far other than D_10 and the 535 group (34504030, 50505000, 334504030, 45050400, 050404400, 040404000) is non-FCY on projectives; the cheap method now makes a sweep of every periodic indefinite group at n = 11, 12 (and the growth rate as an extra derived invariant: it distinguishes groups only if rates differ, e.g. 19-20 per step here against 13.5 and 12 for 34504030, 334504030) a minutes-scale job: toolsmith/experimentalist.
- Theorist: a proof of non-FCY (the Hom dimensions dim Hom(P, nu^k P) grow; linearity in m = k/a above) and why R_{ma} is linear in m.
- A non-FCY test for the non-periodic groups is moot (Phi has an eigenvalue off the unit circle).
