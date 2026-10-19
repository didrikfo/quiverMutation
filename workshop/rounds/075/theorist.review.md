# Review of workshop/rounds/075/theorist.md

referee: scholar · round: 075
verdict: minor revision

## Reproduction

Hand proof; no script. I re-read Sect. 7 (l. 895-1043) and Def. 3.2 (l. 383-400) of the v1 .tex and redid Steps 1, 3, 4 and the sign bookkeeping of Step 5 myself, plus the side example (A_3/rad^2, i = 2). All of it matched the submission.

- hat(p' alpha_*) = 0: with d alpha_* = -alpha-bar, hat(p'alpha_*) = p'alpha_* + (-1)^{|p'|}(-(-1)^{|p'|}) p'alpha_* = 0. Correct. It uses |p| = |p'| since |alpha_*| = 0.
- The relation signs match the paper's printed relation (l. 1014-1021) in all three cases. For t(phi) != i, (-1)^{|phi-bar|} = -(-1)^{|phi|}. For the loop case t(phi) = s(phi) = i, |phi-bar| = |phi| and the extra -1 comes from -\overline{d phi}. Both agree with rho_x via ||phi|| = |phi| - 1 + chi(t phi).
- (sum_b b^* \overline{b phi})/alpha-bar = 0, because the rightmost letter is phi-bar, not alpha-bar.
- The no-overlap claim (Step 2) holds: leading words begin (leftmost) with phi, alpha^* only ever sits to the right of them, and no leading word begins with alpha^*.
- Side example: M = {phi (deg 1), psi phi (deg 0), psi^* (deg 0), omega (deg 1)}, d phi = psi^*(psi phi), d omega = psi phi. H_0 = k[3 -> 2] x k, and phi - psi^* omega is a non-bounding degree-1 cycle. This agrees. A disconnected H_0 cannot be derived equivalent to the connected A_3/rad^2, so "silting, not tilting" is right given Thm 7.5.

## True?

I found no error. Gaps, in order of weight:

1. The step "Theta is an algebra map" (Step 4) rests on hat p hat q = hat(hat p . q) and "hat p . q = pq unless q begins with alpha^*". This is right (the alpha_* tails of hat p meet q's leftmost letter, and alpha_* phi = 0 for phi != alpha^*). But the submission never says that no right factor q of a normal word begins with alpha^*, which is the reason the exception never fires. It is true: alpha^* occurs only as the last factor. One sentence should say it.
2. Step 5, the alpha phi case, silently drops alpha-bar beta^* in d(alpha-bar phi-bar). That term is zero by the relation phi alpha^* = 0 for phi != alpha_*, with phi = alpha-bar. The sign remark covers the other term only. State it.
3. W1 (paper's presentation of Gamma, "Q-bar paths = Q-words", EM = 1) is taken on trust, as flagged. Gamma's presentation is stated in the paper without proof of completeness of the relations. So the "closed" claim in (1) means "closed relative to the printed presentation of Gamma", and the claim should say so.
4. Dimension agreement "in every bidegree" (Step 3) is asserted by a bijection of normal words with paths. I checked degrees for alpha phi (|alpha-bar phi-bar| = ||alpha phi||); the other letter types are not shown.
5. The headline says "dg isomorphism ... (no acyclicity of D(i))". The proof does not use acyclicity anywhere I can see. But the dg-ideal J = (X, dX) being a dg ideal relies on d^2 = 0 (Sect. 8, E-202), which E-202 checks only on 7 082 random letters. The conclusion is therefore conditional on Thm 8.x. The submission mentions this but leaves it out of the title.

## New?

Partly. E-214 lists "Gamma_red = N" as open and unchecked, and E-189/E-199/E-211 cover only the downstream consequences. Grepped `Gamma_red`, `hat(p`, "Gamma_red" in FINDINGS, HYPOTHESES, RETRACTIONS and EXPERIMENTS: it appears only in E-214. Nothing in RETRACTIONS. Steps 3-5 are new to the record. The side finding (silting, non-tilting at i = 2 of A_3/rad^2) is an illustration of E-189/E-151. That is not new, as the author says. The paper's own example (l. 403-436) is the same quiver, but there it is run at a different vertex.

## Evidenced?

Adequate for a proof sketch, weak as evidence. The only worked example has X empty (N = kM = Gamma_red), so it cannot distinguish the claim from a vacuous identity. The step the claim is about (X non-empty, cancelling phi alpha^{-1}) was checked on no example. The author concedes this and proposes a machine check. There is no stated range (arrow degrees, loops at i, parallel arrows). Parallel arrows i -> j are explicitly unchecked (W2).

## Scope

Title: "holds as a dg isomorphism under Thm 1.1's hypotheses alone". Narrowed wording: "Theta: N -> Gamma_red is a bijective dg-algebra map, assuming the paper's presentation of Gamma (l. 905-915), its Lemma 7.4 / Constr. 7.3, and d^2 = 0 on M (Thm 8); checked by hand in one proof sketch, with no example in which X is non-empty." The claim "(2)... closed at the same level of rigour as B of E-214" is fair only with that qualifier.

## Required for acceptance

1. Add the one-line reason in Step 4 why no right factor of a normal word begins with alpha^* (so the "unless q begins with alpha^*" clause is never used).
2. Add the vanishing alpha-bar beta^* = 0 to the alpha phi case in Step 5.
3. Narrow the title and claim (1) to "relative to the printed presentation of Gamma and d^2 = 0 on M", per Scope.
4. Run Theta on one hand example with X non-empty: i = 1 in A_3/rad^2 (the E-214 example). Compute N and Gamma_red in degrees 0 and 1 and report the dimensions side by side. This is doable in one sitting and is the only test of the claim's content. The machine check on random parents is `[next round]`.
5. State the degree check of the bijection in Step 3 for all four letter types (phi, alpha phi, alpha^*, and the t(phi) = i variants), not only alpha phi.
