# Theorist notebook (rewritten round 086)

## What I believe now
- 086 (T5): Lemma A (monomial, X1 and Y disjoint, D = 0 on arrows) has a direct proof: bad sites (theta_xy or adjacent x y) never overlap because X1 cap Y = empty; D preserves the tile skeleton and creates no new site; skeletons with a site are P(x)[theta -> xy](x)Q, contractible; skeleton with no site = normal words of T/(xy). Works over R = k^n. Machine: 300 random quivers n = 1..5, 13 748 cells, 0 mismatches; control with a shared arrow: 164 of 215 fail. So E-235's step 4 no longer cites Anick; theorem number for the original is moot (UNVERIFIED stays for Anick/Halperin-Lemaire).
- 083 (T5): Prop. 7.6 cancellation (kM -> N quasi-iso) proved by hand under non-negative grading + no loop at i of ANY degree. Remaining unproved inputs: "d lowers ||.|| by 1" (machine-checked 15 966 letters, not line by line), v1 LaTeX only. Reduction: omega = ||w|| + #stars, gr d = d_1, then W-count filtration to Lemma A.
- Every `noCycles` parent satisfies the hypothesis. E-211 satisfies it: cancellation fine, failure is step 7.
- OPEN: loops of degree >= 1 at i (d_1 not monomial, 6968 of 10206 loop letters differ).
- Older: G2/step 7 (065, 067): repo kernel = Oppermann's iff D(i) acyclic (E-199). 075: Gamma_red = N (E-217). T7: LNA single-mutation rule proved (E-187, E-205); H-010 supported unproved.

## What I tried
- 086: theorist_lemmaA.py (random graded quivers, exact ranks vs normal-word counts, overlap control). First long run ran into the 10-min limit by weight 7 with 400 instances: size with small runs first.
- 083: theorist_starfilt.py, theorist_omega.py. Earlier: twoby, cycD, Gamma_red, 073, 070, 067, 065, 062.

## Next
- Loops at i of degree >= 1: loop-count filtration or matrix form (psi, theta1, c, theta2).
- Apply the lemma literally to gr_G of the real kM on a parent with n <= 5 (needs signed products) rather than to abstract quivers.
- Derive "d lowers ||.|| by 1" line by line from Def. 3.2; d^2 = 0 beyond E-202; W1/W2, Thm 7.5.
- T7: two pair-sharing bystanders, k = 4 needs OVERNIGHT.
- Blind spots: proofs are relative to the printed presentation; the lemma check is for abstract quivers, not for the E-235 complex; capped enumerations are not laws.
