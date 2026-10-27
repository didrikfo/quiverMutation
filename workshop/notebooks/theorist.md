# Theorist notebook (rewritten round 090)

## What I believe now
- 090 (T10/T5): loops of the mutated quiver M (Oppermann Steps 1-4): (a) loop of Q, (b) alpha.phi from a 2-cycle, (c) phi.alpha^-1 with phi: i->v of degree >= 1 parallel to alpha in A, (d) Step 4. M has a loop at i iff Q has one at i (any Q). My 089 claim "ungraded parent -> no positive-degree loop" is TRUE ONLY for parents with all arrows in degree 0 and no relation arrows. A parent kQ0/(R) enters with a degree-1 arrow omega per relation (const.usual_algebra l. 291); for acyclic such parents the loops of M are exactly the (c) pairs (omega out of i, alpha in A, same target): at v != i, degree >= 1. Referee example 0->1->2, 0->2, omega 0->2: loop at 2 after mutating at 0.
- The 25 failing n = 7 parents (16 class 1, 9 class 2) have acyclic Gabriel quiver, 78 relations, none a cycle (theorist_relarrows.py parents, needs pickles from rounds/050/toolsmith_collect.py, 442 s / 335 s). So no loop at i of any degree: E-235 applies, the 25 failures are not cancellation failures. Positive-degree loop at i only arises if dg mutation is iterated without H_0 (depth 2 with relation arrows, 27 of 1 706 nodes' legal steps; depth 3 hereditary, 089).
- E-216's refusal is on the Gabriel quiver of H_0, not on the dg quiver.
- 086 (T5): Lemma A (monomial, X1 and Y disjoint) proved; 13 748 cells. 083: Prop. 7.6 cancellation proved by hand if non-negative grading + no loop at i; unproved inputs: "d lowers ||.|| by 1" (machine-checked only).
- Older: repo step-7 kernel = Oppermann's iff D(i) acyclic (E-199); Gamma_red = N (E-217); LNA single-mutation rule proved (E-187, E-205).

## What I tried
- 090: relation-arrow Parts A/B (9 020 steps, 0 mismatch with (c) prediction); 25-parents check. 089: one-step classification, iterated scan. 086 lemmaA; 083 starfilt, omega; earlier twoby, cycD, Gamma_red, 073, 070, 067, 065, 062.
- Lesson: I generated test parents in the form I assumed, not in the form the source prescribes. Read the construction of the input before sampling it.

## Next
- Derive "d lowers ||.|| by 1" line by line from Def. 3.2; d^2 = 0 beyond E-202; Thm 7.5.
- Check whether the (c) loops at v != i are essential in a minimal model, and that H_0(kM) = End(T) is unaffected.
- Apply Lemma A to gr_G of the real kM on a parent with n <= 5 (signed products).
- T7: two pair-sharing bystanders, k = 4 needs OVERNIGHT.
- Blind spots: proofs relative to the printed v1; Ext-arrows of acyclic algebras being non-loops used without source check; capped enumerations are not laws.
