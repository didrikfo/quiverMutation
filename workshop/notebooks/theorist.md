# Theorist notebook (rewritten round 075)

## What I believe now
- G2/step 7 (065, 067): with TRUE paths, (A) g = gd always, (B) red(B_1 e_i) = rbar.I and K = {z_0: (d z_0)/a in I all a} always (math). The repo's kernel = (B) iff every tail p is vertex-simple, in particular when D(i) is acyclic (E-199).
- 073: cycle meeting D(i) is NOT harmless (rounds/073/theorist_cycD.py; monomial n = 4, i = 1: repo drops rbar.b.c). Guard inert on cyclic input; `noCycles` (search.py l. 377/420) keeps search parents acyclic, so only direct callers of mutateAtVertex are exposed.
- 075: Gamma_red = N := kM/(X, dX) as dg algebras, under Thm 1.1's hypotheses only (no acyclicity needed): diamond lemma on T(V)/(rho) (no overlaps), basis match with hat p, hat(p' alpha_*) = 0, d matches on generators (sign of alpha phi checked: |alpha-bar| = -1). Weak: paper's presentation of Gamma, hat p hat q = hat(hat p q) taken as printed. So End(T) = H_0(kM) is now: Thm 7.5 + this + E-214 B. By hand only; after review (075 response) stated relative to the printed Gamma presentation and d^2 = 0 on M, with an X-nonempty hand example (A_3/rad^2, i = 1: dim 6 = 6, relation killed, d agrees) and a dimension-only script on 8 cases.
- J = 0 is NOT implied by g = gd (lives in H_1(M)). Premise rests on per-step J_i = 0 check, Thm 1.1 (E-189), left/right convention (E-195), generation (E-168).
- Hand: A_3/rad^2 at i = 2 is silting not tilting (H_1(kM) != 0, H_0 = A_2 x A_1).
- T7: right (062) and left (070) mutation of an LNA row proved relative to the repo's procedure, n = 4..10. Open: two bystanders, k>=4. H-010 supported, unproved; T1, T2 closed; End(T) check blind to J (E-174, E-181).

## What I tried
- 075: Gamma_red = N (hand proof, no script). 073: cycD. 070: leftcheck. 067: cycoff, cycup. 065: g2exact, g2cyclic. 062: t7check. 059: enddim. 054: gen, loops.

## Next
- Machine check of Theta/H_0 on tiny algebras with X non-empty (explicit complexes of projectives vs H_0(kM)); the first independent test of End(T) = H_0(kM).
- Skeptic review of W1/W2 of 075 against the .tex.
- Cyclic D(i) (E-211): decide H_1(kM) for the corrected child (gd 4, four vertices).
- True-path step 7 sweep on random finite-dim cyclic algebras; T7 two-bystander window.
- Blind spots: proofs are about the repo's procedure or a printed presentation, not tilting complexes; one example is not a rate; Coxeter match is necessary not sufficient.
