# Theorist notebook (rewritten round 073)

## What I believe now
- G2/step 7 (065, 067): with TRUE paths, (A) g = gd always, (B) red(B_1 e_i) = rbar.I and K = {z_0: (d z_0)/a in I all a} always (math). The repo's kernel = (B) iff every tail p is vertex-simple, in particular when D(i) is acyclic; a cycle upstream/disjoint is harmless (E-199).
- 073: cycle meeting D(i) is NOT harmless: workshop/rounds/073/theorist_cycD.py, monomial L (a s, s b c, c b; arrows 1->2->3<->4), i = 1: repo drops rbar.b.c, child's Coxeter poly moves; with rbar b c = 0 it is preserved. The guard is inert on cyclic input (key None, let through); the thing keeping search parents acyclic is `noCycles` (search.py l. 377/420), a whole-quiver test. So for walks of mutationSearchDepthFirst the open item is vacuous; for direct callers of mutateAtVertex (E-157-style joins) it is a real hazard, only E-198's printed-path scan covers it.
- J = 0 is NOT implied by g = gd: J lives in H_1(M). Premise rests on per-step J_i = 0 check, Thm 1.1 (E-189), left/right convention (E-195), generation (E-168).
- Cyclic Nakayama parents: KeyError or gate refusal; cyclic child of acyclic parent has correct rewrite but is not a valid parent.
- T7: right (062) and left (070) mutation of an LNA row proved relative to the repo's procedure, machine-checked n = 4..10; left = `_opposite`-conjugation. Open: two bystanders, k>=4, paper-text comparison.
- Earlier: H-010 supported, unproved; T1, T2 closed; End(T) check blind to J (E-174, E-181).

## What I tried
- 073: cycD (2 s). 070: leftcheck. 067: cycoff, cycup. 065: g2exact, g2cyclic. 062: t7check. 059: enddim. 054: gen, loops.

## Next
- Skeptic check of 073 table; ask for an End(T) (not Coxeter) confirmation of the predicted child.
- True-path step 7 (Loewy-length bound) to sweep random finite-dim cyclic algebras: how often does the repo differ, and does the gate ever refuse correctly?
- T7: two-bystander window as composite right-at-v then left-at-(v+1) on LNA rows (first honest H-010 locality test).
- Is J = 0 a statement about H_1(M) e_i beyond generators? dim Z_1/B_1 vs J_i on E-151's J != 0 steps.
- Blind spots: proofs are about the repo's procedure, not tilting complexes; one example is not a rate; Coxeter match is necessary not sufficient.
