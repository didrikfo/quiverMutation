# Theorist notebook (rewritten round 065)

## What I believe now
- G2 (065, T5): for parents with no oriented cycle (parallel arrows fine), repo step 7's kernel EQUALS Oppermann's red(Z_1 e_i) as subspaces of {(r, tail)} (proof: every nonempty q starts with an arrow of A_i, so red keeps the q-empty part z_0; z is a cycle iff (d z_0)/alpha in I for all alpha; Lemma L: red(B_1 e_i) = rbar.I). Same sign, no transpose in the formulas (steps 4, 6 all +, step 5 relative -). Tested 7 479 pairs n = 5..7, 0 differences (theorist_g2exact.py). Counts: repo r7 >= g (= #degree-2 arrows i->v), excess = Lemma-L elements `reduce` drops; g = gd observed, unproved.
- J = 0 is a per-step finite check on the parent (kernel of x -> (x alpha); AI 2.32(b), E-099/E-171 already record it). With Prop 1 the premise reduces to: (a) every parent on a printed path acyclic (UNCHECKED scan), (b) Thm 1.1 cancellation + left/right convention, (c) generation (E-168). More tests cannot close (b); (a) is a scan.
- Cycles: the repo's allPathsBetween is vertex-simple, so cartanMatrix, isInIdeal, _minimalGenerators are silently wrong on cyclic quivers; a cyclic CHILD of an acyclic parent has a correct rewrite (true-path Cartan = R C R^T) but is not a valid parent. Cyclic Nakayama parents: KeyError (_spliceOutVertex splices only the first passage) or gate refusal.
- T7 (062): right mutation of an LNA row passes the gate iff v <= n-1 and no relation ends at v+1; result LNA iff additionally v = n-1 or r_v = 2; Catalan factorisation; relative to the repo's procedure, not the paper (E-187). Left mutation and two bystanders open.
- Earlier: H-010 supported, unproved; T1 closed; T2 closed; T10/T5 End(T) check blind to J (E-174, E-181).

## What I tried
- 065: theorist_g2exact.py, theorist_g2cyclic.py. 062: t7check. 059: enddim, enddim2. 057 kfit. 054 gen, loops. 051 rulelen, ablate.

## Next
- Prove g = gd, or find x with red(z_x) = 0 (would give H_1(M) != 0, i.e. a J-type class in the dg model).
- Extend Lemma L to parents with a cycle away from i; decide whether Prop 1 needs vertex-simple tails or true paths (cyclic parent with a running rewrite not yet found).
- Ask toolsmith for the cycle scan on printed paths; scholar for the cancellation step; skeptic for the APR hand check.
- T7 two mutations: after a slide (r_{v-1}, r_v) = (2, 0), which v' give LNA again (transition graph on LNA rows; attack H-010).
- Blind spots: Step 0's list of M's degree-1 arrows is read from v1 Def. 3.2, not derived; repo-procedure proofs, not tilting-complex statements; tests are random small multigraphs, not walk steps.
