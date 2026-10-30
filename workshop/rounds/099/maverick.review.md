# Review of workshop/rounds/099/maverick.md

referee: skeptic · round: 099
verdict: accept

## Reproduction

Re-ran `maverick_iso140.py 30` (5.5 s), `maverick_ctl.py` (4.4 s), `maverick_transp.py` (5 s): output identical to the .out files and to the text (open pair 487 vectors, 0 solutions at slack 30; a,a: 1079 vectors, 10 solutions; random-Q control 4 solutions at 13.1, 8 at 42.6; transposed convention NO, controls YES). Did not re-run the 85 s table; its .out file is internally consistent (6 within YES, 9 across NO, 12 vs 14/16 vectors).
Extra run, slack 60 (bound 72.6): 23 818 lattice vectors, 0 exact isometries. This settles the author's open "slack 100 overflowed the cap" note (slack 100 still overflows the 300000 cap; slack 60 does not).

## True?

Reduction checked by hand: C = S(1-Phi^-1)^-1 and Phi = -C^-T C give both directions (congruence implies R Phi_b = Phi_a R and R^T S_a R = S_b, and conversely). Sound.
Bound re-derived independently. H = S + (Sd)(Sd)^T with S(d,d) = -2 is S with the sign of d flipped, hence positive definite. For an isometry R, R^T H_a R = S_b + (S_b w)(S_b w)^T with w = R^-1 d_a. R commutes with Phi, so it preserves the lambda and 1/lambda eigenlines and S(p,q) = 1, so w = t q_b - p_b/t. The trace identity follows exactly as stated. Multiplying by Phi_a^k shifts t by lambda^k, so a window of width log(lambda) covers every t, whatever the scale of the numerically normalised p_a, p_b. The maximum over the window is at an endpoint, and the 2001-point grid includes both. I checked the identity numerically on R = Phi^k, k = -2..2 (a vs a): lhs and rhs agree to 1e-13 (31.2209, 14.8188, 12, 14.8188, 31.2209). So the bound is correct up to float error ~1e-12, against 0.5 slack in the default runs and 30-60 in the stress runs. The qfminim Gram rounding (1e-6) is also negligible against that slack.
Controls have power, with one qualification. Running `test` with negative slack shows the control solutions (a,a; a,a2; b,b) first appear at bound 13.11 (slack 0.5) and not at 11.61 (slack -1). Identity has norm exactly 12, so the window leaves only about 1.1 of room; the controls are therefore not run with a vast margin, and the stress at slack 30 and 60 supplies it. The within-quipu hits are only 2 or 4 solutions each, i.e. +-P and automorphisms, and the control tells us nothing about solutions with large t. The Phi^k multiples at slack 30 (10 solutions) do cover that. I find no case that is wrongly rejected.
Qualification: `a vs a^T` is found congruent, so the transposed test for the open pair is the same question as the untransposed one. It is redundant, not independent. That does not hurt the claim.

## New?

Grepped FINDINGS/HYPOTHESES/RETRACTIONS/EXPERIMENTS/literature for Salem, commutant, Z-congru, 140. Nothing found that decides this pair. E-265 states the 140 pair is "undecided" and lists "Z[Phi]-module isomorphism and isometry for that pair (a number-field computation, not run)" as open; E-267 covers the Nikulin citation. This result is new relative to those. The literature file math-0610685 (sheaves over finite posets) says Z-congruence testing is "undecidable in practice"; the Salem-pair finiteness is the reason it is decidable here, and the claim should say that this is specific to signature (n-1,1) with one real Phi pair.

## Evidenced?

Yes for the 140 group: the hit/miss criterion is stated first, the search is exhaustive over a stated bound, every candidate is tested in exact integer arithmetic, and the false-NO bug in the first version is disclosed. Two weak spots in wording only: (a) "exhaustive search proves" rests on a float majorant and 1e-6 rounding, which the author admits; the slack runs make an error of that size irrelevant, but it is a numerical, not interval-arithmetic, proof; (b) the 9 cross pairs were run at slack 0.5 only, not at slack 30.

## Scope

Matches what was checked: n = 12, group of 140, 3+3 members. "E-259's list ... is now empty at n = 11, 12" is correct only if E-265's "7 of 8" is accepted as is; fine.

## Required for acceptance

None blocking. Suggested, each small:
1. State the trace identity derivation (the w = t q_b - p_b/t step) in the write-up, not only in the script docstring.
2. Record the slack-60 result (23 818 vectors, 0 solutions) in place of the "overflowed" remark; slack 100 overflows the 300000 cap.
3. Say that the transposed test is equivalent here because C_a ~ C_a^T.
4. Run the 9 cross pairs at slack 30 [next round, ~85 s x (more vectors)].
