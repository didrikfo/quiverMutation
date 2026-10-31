# Theorist notebook (rewritten round 103)

## What I believe now
- 103 (T5, literature): a loop of degree >= 1 at i in a Construction-2.2 dg model forces a closed Q0-path through i (induction on degree: d phi is a nonzero cycle of smaller arrows), so i is on an oriented cycle and D(i) is cyclic. `mutateAtVertex` refuses that and the search expands only acyclic nodes (E-223). So Oppermann's Step 4 is empty for every library step and E-235's "no loop at i" holds automatically there. The converse (cycle through i => loop) is not claimed.
- Hypotheses matched: Oppermann Thm 1.1/7.5/7.6 need only no degree-0 loop (l. 263, 316, 322; reading: degree-0 loops would put P_i into the approximation); AI Thm 2.31/2.32 and Keller math/0601185v5 Thm 3.12 (Rickard) have no quiver hypothesis (numbers by counter order, Rickard's own number UNVERIFIED). arXiv has only Oppermann v1.
- What fails for a loop phi of degree >= 1 is Lemma A's "D theta = xy": d_1(phi alpha^-1) = +-phi alpha* + sum beta*(beta phi alpha^-1) + lower. Our gap, not the paper's. Plan (untested): filter by the number of Step-3 letters, kill Step-4 letters on gr, then Step-3 letters on E^1.
- Library on cyclic inputs: 2 x 300 random monomial cyclic algebras: gate True with i on a cycle in 495 / 491 cases, default refuses all, `allowCyclic=True` is KeyError in 394 / 401 and returns (unproved, vertex-simple paths) in 101 / 90. No Hom check of those children.
- Older: 101 cancellation is Prop. 7.6, Anick numbers secondary (Gartner). 098 (T4) class A/B decided by first step (E-187/E-205). 093 conditional on omega-reduction, d_1 in Lemma A (600 parents, 0 failures).

## What I tried
- 103: read Oppermann l. 252-400 and 895-1043, AI l. 340-1040, Keller ICM Sect. 3; two scripts (theorist_loopdeg.py, theorist_cycsweep.py). Lesson: check which vertex a relation sits on (c b is a loop at 4, not 3) and read the path-composition order before building a parent.

## Next
- T5: two-stage elimination on one explicit parent with a degree-1 loop (needs Ext^2(S_i,S_i) != 0) and an independent End(T); Hom-verify one `allowCyclic=True` child with i on a cycle (skeptic).
- T4: shape of the 14 width-5 middle algebras; the 6 step-2 failures (w6/w7); one-step rule for the width-5 class-B rules (E-266).
- Blind spots: Construction 2.2 minimality assumed in the loop argument; Rickard 1989 and Anick unread at source; sweep is capped and random.
