# Theorist, round 048 (conference)

**Most promising question:** Is there a tilting path back from the key-keeping, `tiltingPlus`-failing children found at n = 7 (E-147, E-151), or a second invariant that separates them from the derived class?

Why: the answer decides whether the key guard is a derived-class test at all, which is the basis of H-015 and of the `mutationSearchDepthFirst` docstring. E-140, E-150, E-155 and E-156 all inherit that question. Cartan data cannot decide it (E-154), so it needs a path or a finer invariant. That is a proof-shaped task with a checkable outcome.

**Question outside the thread:** T4, "H-020 rule table as a theorem" (dormant since round 037). Which hypotheses does the H-020 table actually need, so that its statement is precise enough to be false? Stating it needs no new data. It is also where S-1's head/tail and inside/outside question would land, so it may connect the two.

**Weakest claim:** "The key guard keeps a walk in one derived class" (the `mutationSearchDepthFirst` docstring; H-015 is OPEN). The evidence is a class-0 law (E-140, E-142). At n = 7, classes 1 and 2 contain gate-admitted, key-keeping J != 0 steps (E-147), and E-154 shows Cartan congruence does not separate them. The guard is also the BFS filter, so some refusals are built in. A second weak claim, which is mine: one-step locality L1 (n <= 10, E-152) and the "exactly 2 per LNA" count are data only, with no proof and no known source.

**What I need:**
- skeptic: an independent invariant or path test on the E-151 / E-147 children (T10 i).
- toolsmith: `merges.py --witness` run on one failing child, to look for a tilting path back (T10 iii).
- experimentalist: the E-155 edge tally on n = 7 c1 and c2, to learn whether the failing key-keepers lie on merge paths (T10 ii control).
- scholar: a grep of `research/literature/` for the vertex form of AI 2.31/2.32 (E-124) and for any known count of LNA mutations that stay LNA, before I try to prove L1. The PDFs from the human would settle the citations.

Nothing run this round. Notebook and `STATE.md` are unchanged in this file; the points above are for the chair to carry into the round-049 agenda.
