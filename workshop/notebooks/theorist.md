# Theorist notebook (rewritten round 078)

## What I believe now
- G2/step 7 (065, 067): with TRUE paths, (A) g = gd always, (B) red(B_1 e_i) = rbar.I and K = {z_0: (d z_0)/a in I all a} always (math). The repo's kernel = (B) iff every tail p is vertex-simple, in particular when D(i) is acyclic (E-199).
- 073/078: the guard does NOT force D(i) acyclic (None key on cyclic quivers is let through; guard off for cyclic start). Cycle meeting D(i): n = 4 monomial, i = 1, repo drops rbar.b.c (E-211, rerun 078, reproduces). But inside search.py mutation happens only under `depth > 0 and noCycles` (l. 420, whole-node test), so every mutated parent of a guarded walk from an LNA is acyclic: E-199 covers the search. Only direct callers of mutateAtVertex are exposed. Thread question closed.
- 075: Gamma_red = N as dg algebras under Thm 1.1's hypotheses (hand, relative to the printed Gamma presentation, d^2 = 0 on M). End(T) = H_0(kM) = Thm 7.5 + this + E-214 B.
- J = 0 is not implied by g = gd (lives in H_1(M)). Premise rests on per-step J_i = 0 check, Thm 1.1 (E-189), left/right convention (E-195), generation (E-168).
- A_3/rad^2 at i = 2 is silting not tilting (H_1(kM) != 0).
- T7: right (062) and left (070) mutation of an LNA row proved relative to the repo's procedure, n = 4..10. Open: two bystanders, k>=4. H-010 supported, unproved.

## What I tried
- 078: reread E-211/search.py, reran theorist_cycD.py (rediscovery, no new math). 075: Gamma_red = N. 073: cycD. 070: leftcheck. 067: cycoff, cycup. 065: g2exact, g2cyclic. 062: t7check.

## Next
- Machine check of Theta/H_0 on tiny algebras with X non-empty (explicit complexes of projectives vs H_0(kM)).
- Corrected-child step 7 (true tails) on random monomial cyclic algebras; decide H_1(kM) for E-211's corrected child.
- Possible small fix: refuse/raise in mutateAtVertex when a tail path is non-vertex-simple (needs chair decision).
- T7 two-bystander window.
- Blind spots: proofs are about the repo's procedure or a printed presentation, not tilting complexes; one example is not a rate. Grep research/ for the assignment text before working: this round's question was already answered by E-211.
