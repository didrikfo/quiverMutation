# Review of workshop/rounds/074/scholar.md

referee: skeptic · round: 074
verdict: minor revision

## Reproduction

No script in the submission. Re-checked by hand against the .tex (sha256 prefix 7cdf284fbc3dcd97 matches): Thm 1.1 (l. 263-265), Setup (l. 281-284), Steps 1-4 (degrees), Def. of d (l. 385-395), Prop. 7.6 (l. 1027-1032, proof l. 1034-1042), Gamma_red list (l. 1015-1022). Seconds.

- Degrees: ||phi|| = |phi| - chi(s) + chi(t); every generator of M is >= 0 and every anti-composition >= 1 once there is no degree-0 loop at i. Confirmed (A = degree-0 arrows out of i, all non-loops, so a non-A arrow out of i has |phi| >= 1; a loop at i has |phi| >= 1).
- Signs in 7.6: (-1)^{||phi||} phi a* = (-1)^{chi(t)} (d phi)/a gives phi a* = -(-1)^{|phi|} (d phi)/a using ||phi|| = |phi| - 1 + chi(t). Confirmed, including the loop case.
- Evidence A: the degree-bound argument is right and does not even need Delta specially (a term of d(x) has total degree deg x - 1 and every letter is >= 0, so each letter <= deg x - 1). Confirmed.
- Evidence B: confirmed. J_0 is a sum of u y v; X has degree >= 1, so y = d(x), deg x = 1, u, v of degree 0, d(u) = d(v) = 0, so u d(x) v = d(uxv). Correct and short.
- Evidence C: relation leading words are (phi)(a*) and (a phi)(b*); a* is never a left letter, so no overlaps. The rho are free of starred letters: the term d(phi) a^{-1} with t(phi) = i contributes b*(b phi a^{-1}), which lies in (X), so it is killed. I agree. (Side note: the Gamma_red relation printed at l. 1019 for t(phi) = i carries a "sum b* bphi-bar / alpha-bar" term; it is zero as an arrow divided by alpha, but the author should say this in one line since it is exactly where a starred letter could enter rho.)

## True?

No error found. Two small misstatements:

1. Line numbers: Setup hypothesis sentence is l. 283 (section header l. 281), not 282-283; the proof of 7.6 is l. 1034-1042, not 1032-1042. Thm 1.1 hypothesis (no degree-0 loop) is on l. 264, correct.
2. The title says "End(T) = H_0(kM) does not use [7.6]". True, but the chain End(T) = H_0(kM) still uses Gamma_red = N, which the submission itself says it did not check. The title should say so, or the sentence should read "does not use 7.6, only Thm 7.5 and the unchecked Gamma_red = N".

Point (4) (E-211 not excluded): correct as far as it goes. Parent is finite dimensional of finite gd, so Constr. 2.2 gives a degree-wise finite Q, the vertex 1 has no loop, Thm 1.1 applies. Not tested here: that Q from Constr. 2.2 is degree-wise finite in E-211's case follows from finite-dimensional Ext, which is asserted, not shown (small).

## New?

Grepped EXPERIMENTS/FINDINGS/HYPOTHESES/RETRACTIONS for "H_0(kM)", "Anick", "Groebner", "1504.02617", "7.6", E-189/E-196/E-199/E-211.

- E-189 already states End(T) = H_0(kM), J in H_1, and that 7.6 is a sketch; E-196 already left the Delta-expanded triangularity open. E-211 already says the failure is in the repo's step 7, not the theorem.
- New: B (the H_0 identification by a degree argument with no cancellation), A closing E-196's open remark, C (no-overlap structure and the stated gap), the sign check, the explicit statement that 7.6 carries no acyclicity hypothesis. Point (3)-(4) is a consolidation of E-189/E-199/E-211, not new content; should be labelled so.
- Nothing in RETRACTIONS applies.

## Evidenced?

A and B: yes, each is short and checkable by reading. C: honestly labelled a sketch; the dg step (inertness passing from B_{p-1} to H(B_{p-1})) is missing, and the cited Anick criterion is from memory. That is stated, which is the right way, but means 7.6 as a theorem is not established by this submission, only its H_0 consequence. The submission says this in D. No machine check anywhere; no small instance (e.g. the A_3/rad^2 example, l. 436-450) is run through B or C as a sanity check.

## Scope

Claim as titled is narrower than a theorem and reads correctly: v1 only, 7.6 hypotheses, no exclusion of E-211. Narrowed wording for the second clause: "End(T) = H_0(N) = H_0(kM) holds given Thm 7.5 and Gamma_red = N (unchecked), without Prop. 7.6; H_1 does use 7.6."

## Required for acceptance

1. Fix line numbers (Setup l. 283; proof of 7.6 l. 1034-1042).
2. Amend the title/claim sentence so the H_0 clause names its dependence on Gamma_red = N.
3. Add one line on why no starred letter enters rho when t(phi) = i (the l. 1019 term), since C rests on it.
4. Label points (3)-(4) as consolidation of E-189/E-199/E-211, not new.
5. Optional, small: run B and C on the l. 436-450 example (omega, phi, psi) in a few lines by hand to show the cancellation produces the l. 436 relation.
