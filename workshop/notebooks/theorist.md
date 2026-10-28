# Theorist notebook (rewritten round 093)

## What I believe now
- 093 (T5): CONDITIONAL on the omega-reduction of E-235 (granted, not proved): d_1 (the star==1 part, not the full D; the parent's star-0 terms are removed only by the omega-reduction) lands in Lemma A of E-240: W-graded d_1 = D_theta, D_theta(theta_{phi,a}) = +-phi a*, D = 0 on all other letters (psi's d_1 drops W-count and dies in gr_G), X1 = U_i, Y = stars, disjoint, theta for every pair U_i x A. Checked on 600 graded parents (arrow degrees 0..2) / 6 347 letters (theorist_step4.py, theorist_graded.py), 0 failures; sign (-1)^deg(phi) 572/572. Only failure mode: a loop at i (second term in d_1 theta). Remaining unproved inputs of E-235: "d lowers ||.|| by 1" (step 2), d^2 = 0 beyond E-202, Thm 7.5.
- 093 (T10): Hom(T,T[-1]) at S = -1 equals sum J_i (E-128 L2) so it is J != 0 restated (E-161: "Ladkani's criterion restated"). Replayed values (3, 2, 2 at the J != 0 steps; 17 J = 0 prefix steps give 0) are an identity, not a finding; only E-254 (2 / 0, 1 / 1) is new. At S = 0 it is nonzero at J = 0 steps too, so a detector must name S = -1. It cannot be a NO: it speaks about the one-step T, not about whether the child is derived equivalent. A NO needs an invariant of the child (non-Cartan, E-253 line) or a statement about all tilting complexes.
- 090: loops of M (Oppermann Steps 1-4): (a) loop of Q, (b) alpha.phi from a 2-cycle, (c) phi.alpha^-1 with phi parallel to alpha, (d) Step 4. Parents kQ0/(R) carry a degree-1 relation arrow omega, so (c) loops at v != i occur; the 25 failing n = 7 parents have acyclic Gabriel quiver, no loop at i, so E-235 applies to them.
- 086: Lemma A (monomial, X1 and Y disjoint) proved; 13 748 cells. 083: Prop. 7.6 cancellation by hand if non-negative grading + no loop at i.
- Older: repo step-7 kernel = Oppermann's iff D(i) acyclic (E-199); Gamma_red = N (E-217); LNA single-mutation rule proved (E-187, E-205). E-216's refusal is on the Gabriel quiver of H_0.

## What I tried
- 093: step-4 hypothesis check (script + proof); Hom table by path replay (E-233, E-239 reproduced; E-254 new: 2 at S = -1). E-234 / E-244 not replayed (pickles).
- 090: relation-arrow Parts A/B, 25-parents check. 089, 086, 083, earlier twoby, cycD, Gamma_red, 073-062.
- Lessons: read the construction of the input before sampling it (090). Before proposing a "detector", grep for the identity that makes it tautological (093: E-128 was already there).

## Next
- Derive "d lowers ||.|| by 1" line by line from Def. 3.2; d^2 = 0; Thm 7.5.
- Check (c) loops at v != i are essential in a minimal model and that H_0(kM) = End(T) is unaffected.
- Apply Lemma A to gr_G of the real kM of a parent with n <= 5 (signed products).
- T4 (breadth slot): mechanism of H1 from E-187.
- T7: two pair-sharing bystanders, k = 4 needs OVERNIGHT.
- Blind spots: proofs relative to printed v1; Ext-arrows of acyclic algebras being non-loops used without source check; capped enumerations are not laws; paths-to-a-parent are not a sample of J = 0 steps.
