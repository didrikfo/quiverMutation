# Theorist notebook (rewritten round 081)

## What I believe now
- G2/step 7 (065, 067): with TRUE paths, (A) g = gd always, (B) red(B_1 e_i) = rbar.I and K = {z_0: (d z_0)/a in I all a} always (math). Repo's kernel = (B) iff every tail p is vertex-simple, e.g. D(i) acyclic (E-199).
- 073/078: guard does not force D(i) acyclic; E-211 is the cyclic counterexample for direct mutateAtVertex callers; search.py only mutates acyclic parents, so E-199 covers the search. Closed.
- 075: Gamma_red = N as dg algebras under Thm 1.1 (hand, relative to printed presentation). J = 0 not implied by g = gd. A_3/rad^2 at i = 2 is silting not tilting.
- T7: single-mutation LNA rule proved for right (062, E-187) and left (070, E-205) relative to the repo's procedure, n = 4..10. H-010 SUPPORTED, unproved.
- 081 (T7, two bystanders, k = 3, margin 2, pair (8:3)(9:3)): 16/16 admissible cells with both bystanders sharing <= 1 arrow (m <= 3, t 6..12) stay at overlap 2, incl. 3 chained cells; 4/4 controls with a run of three plus a second bystander reach 0. Weak cells: (6,3)(11,3) and (6,3)(12,3) reach 0 and 1 LNAs.

## What I tried
- 081: theorist_twoby.py (+_rest, _ctrl). ~100 s per k = 3 cell, which kills larger sweeps.
- 078: cycD rerun; 075 Gamma_red; 073; 070 left check; 067; 065; 062.

## Next
- Two pair-sharing bystanders on opposite sides ((7,3)+(10,3)), and m = 4, t outside 6..12.
- k = 4 needs a faster localMutationSequences or an OVERNIGHT job; propose ~20 cells.
- Find why a run of three lowers (mechanism), via the E-187 gate rule applied to a 3-step sequence.
- Machine check of Theta/H_0 on tiny algebras; corrected-child step 7 on random cyclic monomial algebras.
- Blind spots: proofs are relative to the repo's procedure; capped enumerations are not laws. Grep research/ for the assignment first.
