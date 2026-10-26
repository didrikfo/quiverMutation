# Theorist notebook (rewritten round 089)

## What I believe now
- 089 (T10/T5): classification of loops of the mutated quiver M (Oppermann Steps 1-4): a loop of M at v is (a) a loop of Q, (b) alpha.phi for a 2-cycle i-v (degree |phi| - [v=i]), (c) phi.alpha^-1 for parallel i->v arrows with phi of degree >= 1, (d) Step 4 on a positive-degree loop at i. So from an ungraded parent M has only degree-0 loops (2-cycles), acyclic parent -> loopless. The repo passes to H_0 each step, so its parents are ungraded: "loops of degree >= 1 at i" never occurs for the 25 failing children (submission 089). Iterated dg mutation (no H_0) does create degree-1 loops at depth 3 (A_3 1<-0->2, sequence 1,0,2,0; quiver level only, minimality unchecked). Script theorist_loopreach.py.
- 086 (T5): Lemma A (monomial, X1 and Y disjoint, D = 0 on arrows) proved directly; 13 748 cells, 0 mismatches; E-235 step 4 no longer cites Anick.
- 083 (T5): Prop. 7.6 cancellation proved by hand under non-negative grading + no loop at i of any degree. Unproved inputs: "d lowers ||.|| by 1" (machine-checked, not line by line), v1 LaTeX only.
- Older: repo step-7 kernel = Oppermann's iff D(i) acyclic (E-199); Gamma_red = N (E-217); LNA single-mutation rule proved (E-187, E-205).

## What I tried
- 089: one-step loop classification (2 574 random steps, 0 outside the list); iterated dg-mutation scan from acyclic starts, depth <= 4.
- 086: theorist_lemmaA.py. 083: theorist_starfilt.py, theorist_omega.py. Earlier: twoby, cycD, Gamma_red, 073, 070, 067, 065, 062.

## Next
- Derive "d lowers ||.|| by 1" line by line from Def. 3.2; d^2 = 0 beyond E-202; W1/W2, Thm 7.5.
- If anyone iterates dg mutation: track the differential on the A_3 example to see whether the depth-3 loop is essential.
- Apply Lemma A literally to gr_G of the real kM on a parent with n <= 5 (needs signed products).
- T7: two pair-sharing bystanders, k = 4 needs OVERNIGHT.
- Blind spots: proofs are relative to the printed presentation; classification read off Steps 1-4 (degree of phi.alpha^-1 unverified by a second reader); capped enumerations are not laws.
