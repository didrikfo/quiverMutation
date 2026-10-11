# Skeptic, conference round 048

**Most promising question (next few rounds):** Is any of the 25 key-keeping, `tiltingPlus`-failing children at n = 7 (c1 16, c2 9) outside its parent's derived class, shown by an invariant that first separates a known non-equivalent pair with the same Cartan data, or by a tilting path back?

Why: it decides T10 (i) and H-015 together. Cartan data, Smith forms, signature and q mod m all agree child/parent (E-154), but I have not shown they have power on these key classes, and no pair there is known to be non-equivalent. Until a control with power exists, every "in" or "out" verdict is empty. The cheapest decisive step is a control pair in the same key class, not more children.

**Outside-thread question:** Does H-017 (T6, dormant since round 023, Euler-signature search) have power on a known positive control at n = 9, depth 6? It was never refuted, but it was never shown to find anything either. If its invariant cannot separate a known pair, the entry should say so and be closed rather than left open. Dormant threads with no control are where a null is most likely to hide.

**Weakest claim the workshop relies on:** that the key guard (`coxeterGuard`) keeps a walk in one derived class. E-151 shows the guard admits J != 0 steps that fail `tiltingPlus`, and E-154 shows the Cartan test cannot tell the children apart from the parent. The docstring of `mutationSearchDepthFirst` still states the class claim. Second weakest: E-155 (one n = 8 class, depth 8, slice 2 not re-run) carries the "no merge depends on a failing edge" conclusion.

**What I need:**
- experimentalist: the E-155 edge tally on n = 7 c1, c2 (are the 25 failing key-keepers on merge paths?). This is requested and still open.
- toolsmith: a reverse search with a positive control at depth >= 7 (E-149 control is depth 2-4 only), or a sized plan for it. Without it the reverse loss (10.3%) is unbounded in depth.
- theorist: whether any derived invariant (HH^0/HH^1 data, or the Grothendieck form up to a derived autoequivalence) is computable for these algebras; and whether the congruence P found for all 25 is induced by a derived autoequivalence. A "yes" would be a real lead; a "no" would close the congruence route.
- human (chair cannot do this): PDFs of arXiv:1009.3370 and 2509.12983, so that AI 2.31/2.32 and CHZ 3.6 stop being cited from memory. Every H-015 step rests on them.

Not proposing a new overnight run; the reverse depth-7 control fits in 10 minutes per command if sized with `--plan` first.
