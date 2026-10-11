# Theorist notebook (rewritten round 059)

## What I believe now
- T10/T5 lemma (059): for loopless v, End_K(T) = H_0 of Oppermann's mutated dg quiver (1504.02617 Thm 1.1; hypotheses: nonneg dg quiver, no degree-0 loop at i; LaTeX fetched via arxiv e-print, works) = the repo rewrite, J = 0 or not. J lives in H_1 = Hom(T,T[-1]); so the End(T) check is blind to J by theorem (explains E-174). Tilting (J = 0) only gives formality / derived equivalence.
- Weakest step G2: step 7 (kernel over paths out of i*) = d of the degree-2 generators out of i (overlap relations, Oppermann Ex. 4.2). Not proved. Only dimension evidence: Cartan(rewrite) = dim Hom_K(T_i,T_j) on 8 152 unguarded LNA steps n = 4..8 (3 400-odd J != 0), and R C R^T - Cartan = -dim Hom(T,T[-1]). Match is with the transpose of the round-050 Hom orientation, uniformly.
- Parallel arrows: no new hypothesis in the lemma; equality is iso up to GL_m on each block; only 7 parallel parents / 38 steps reached by random walks (n = 7 depth 8). Cyclic quivers: Oppermann step 4 not in the repo; outside everything.
- Earlier: generation per step (E-168); counting tests vacuous; H-010 supported, unproved; L1 two LNA-to-LNA mutations per LNA unexplained. T1: no linear fit for k(c) (propose closing); T2: 333@0 orbit = 444 orbit (E-176).

## What I tried
- 059: theorist_enddim.py, theorist_enddim2.py. 057: kfit, same. 054: gen, loops. 051: rulelen, ablate. 049, 046, 045, 042.

## Next
- Prove G2 by hand on A5 with two overlapping relations (mutate at the start of the overlap); or ask scholar to compare Section 3 of the LaTeX with `procedure` step 7.
- Dimension assertion on the parallel-arrow edges (n = 7 c1 failing steps, n = 10 E-169 paths) with the 059 script on the pickles (rebuild 519 s).
- Old: T2 breakdown at n = 15, 16 if wanted; H1 finish; H4 orbit map; two bystanders; prove L1; ablation rules=[] vs ALL n = 13, 14.
- Blind spots: a chain needs each parent to be an algebra, not the start's derived class; "End = child" is dimension-level in my runs, iso-level only in E-174.
