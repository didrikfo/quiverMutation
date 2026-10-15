# Theorist notebook (rewritten round 067)

## What I believe now
- G2/step 7 (065, 067): with TRUE paths, (A) g = gd always (proof: red kills rad Z_1; a cycle with q nonempty splits by first arrow a, d(a z_a) = a d z_a), (B) red(B_1 e_i) = rbar.I and K = {z_0: (d z_0)/a in I all a} always. Repo = true wherever descendant set D(i) is acyclic; a cycle upstream/disjoint is harmless (7 389 pairs, 0 differences, theorist_cycup.py). Cycle meeting D(i) is the only open part (repo vertex-simple, E-193).
- J = 0 is NOT implied: g = gd only rules out generators with red(z_x) = 0; J lives in H_1(M). Premise rests on: per-step J_i = 0 check (E-099/E-171), Thm 1.1 cancellation (E-189), left/right convention (E-195), generation (E-168), and the scan "no printed path step has a cycle meeting D(i)" (UNDONE, weaker condition than global acyclicity).
- Cyclic Nakayama parents: KeyError (_spliceOutVertex) or gate refusal; cyclic CHILD of acyclic parent has correct rewrite (true-path Cartan = R C R^T) but is not a valid parent.
- T7 (062): right mutation of an LNA row passes gate iff v <= n-1 and no relation ends at v+1; LNA iff also v = n-1 or r_v = 2; Catalan factorisation (E-187), relative to the repo's procedure. Left mutation, two bystanders open.
- Earlier: H-010 supported, unproved; T1, T2 closed; End(T) check blind to J (E-174, E-181).

## What I tried
- 067: theorist_cycoff.py (low power), theorist_cycup.py. 065: g2exact, g2cyclic. 062: t7check. 059: enddim. 054: gen, loops.

## Next
- Cycle in D(i) with finite-dim algebra (i -> cyclic Nakayama): true-path K vs repo; can the repo's raw kernel even be built.
- Is J = 0 a statement about H_1(M) e_i = Z_1/B_1 in all of degree 1, not only the generators? Try: dim Z_1(i,v)/B_1 vs J_i on small examples (the J != 0 steps of E-151).
- Ask toolsmith for the D(i) scan; skeptic to hand-check (A); scholar the cancellation step.
- T7: after a slide (r_{v-1}, r_v) = (2, 0), which v' give LNA again; left mutation by duality.
- Blind spots: Step 0 list of M's degree-1 arrows read from v1 Def. 3.2; proofs are about the repo's procedure, not tilting complexes; random small multigraphs, not walk steps.
