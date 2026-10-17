# Theorist notebook (rewritten round 070)

## What I believe now
- G2/step 7 (065, 067): with TRUE paths, (A) g = gd always, (B) red(B_1 e_i) = rbar.I and K = {z_0: (d z_0)/a in I all a} always. Repo = true wherever descendant set D(i) is acyclic; a cycle upstream/disjoint is harmless (7 389 pairs, 0 differences, theorist_cycup.py). Cycle meeting D(i) is the only open part (repo vertex-simple, E-193).
- J = 0 is NOT implied by g = gd: J lives in H_1(M). Premise rests on: per-step J_i = 0 check (E-099/E-171), Thm 1.1 cancellation (E-189), left/right convention (E-195), generation (E-168), and the UNDONE scan "no printed path step has a cycle meeting D(i)".
- Cyclic Nakayama parents: KeyError (_spliceOutVertex) or gate refusal; cyclic CHILD of acyclic parent has correct rewrite but is not a valid parent.
- T7: right mutation (062, E-187) and LEFT mutation (070) of an LNA row proved relative to the repo's procedure, machine-checked n = 4..10. Right: gate iff v<=n-1 and no relation ends at v+1; LNA iff v=n-1 or r_v=2; (r_{v-1},r_v):=(2,0); count C_{v-1}C_{n-1-v}. Left is the mirror (slide right; gate: v>=2, no relation starts at v-1; count C_{v-2}C_{n-v}). Left in the repo IS `_opposite`-conjugation, so the proof is formal. Open: two bystanders, k>=4, paper text comparison (2112.08129 is the human's own paper; procedure.py = its steps 1-7 per FINDINGS l.2599; LaTeX not in sources/).
- Earlier: H-010 supported, unproved; T1, T2 closed; End(T) check blind to J (E-174, E-181).

## What I tried
- 070: theorist_leftcheck.py (left mutation direct, 2 min). 067: cycoff, cycup. 065: g2exact, g2cyclic. 062: t7check. 059: enddim. 054: gen, loops.

## Next
- T7: two-bystander window as composite right-at-v then left-at-(v+1) on LNA rows: which composites return an LNA, compare with E-152 cells; this is the first honest test of H-010 locality beyond one step.
- Cycle in D(i) with finite-dim algebra (i -> cyclic Nakayama): true-path K vs repo.
- Is J = 0 a statement about H_1(M) e_i = Z_1/B_1 beyond generators? dim Z_1(i,v)/B_1 vs J_i on the J != 0 steps of E-151.
- Ask scholar for the 2112.08129 e-print comparison of steps 3-7; skeptic hand-check of (A).
- Blind spots: proofs are about the repo's procedure, not tilting complexes; random small multigraphs, not walk steps; a "formal by duality" claim is only as good as the coordinates of the mirror (hence the direct run).
