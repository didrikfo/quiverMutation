# Review of workshop/rounds/054/theorist.md

referee: skeptic · round: 054
verdict: minor revision

## Reproduction
Ran `theorist_gen.py lna 5|6|7` (112/420/1584 steps) and `path13` (13 edges plus the 16 failing steps). Counts match the table exactly: all H1-H3 hold, det = +-1, Hom(T,T[1]) = 0 everywhere, Hom(T,T[-1]) = 0 on 70/252/924 and 13/13, nonzero on all 16 failing steps (dims 1 x14, 2 x2). The saved path13 output differs from mine only by the `time` lines. I did not re-run the 519 s collection. The 16 failing steps come from a pickle I did not rebuild.

## True?
The proof is correct. With no loop at v, every arrow v->h has h != v, so each P_h is a summand of T and M is in add(T). The triangle P_v -> M -> T_v puts P_v in thick(T), and the other P_i are summands, so thick(T) contains thick(A) = K^b(proj A). Parallel arrows only give multiplicities. The argument also covers the empty-M case, where T_v = P_v[1] (the claim says "every vertex with an out-arrow", but out-arrows are not needed). The condition is exactly "no loop at v", and nothing more is used.

Gaps:
1. "So for all tiltingPlus steps on acyclic algebras" is only valid if every intermediate algebra on a path is loopless at the mutated vertex. The author concedes the script asserts this for 29 edges only. A mutated quiver can create loops or oriented cycles, although the quiver here is acyclic. Until acyclicity of every child is checked or argued, the claim "for all steps" is conditional. The proof itself is fine.
2. The hypotheses in the table (H1-H3) are exactly the proof's premises, so the script tests nothing the proof does not already give. It only confirms that the code's data satisfy them. The scripted numbers that carry real content are the Hom(T,T[+-1]) columns. This is acknowledged for det but not for H1-H3.
3. Hom(T,T[1]) = 0 "never fails" is empirical (2145 steps) and is not proved. The text should not suggest otherwise. Also, 8 parallel-arrow failing steps are outside the Hom code (the author lists this in Next).

## New?
Partly recorded. `research/literature/rickard-morita-theory-derived-categories.md` (the generation paragraph, about line 78) already says that for T = mu^-_{P_i}(A), generation holds "because the triangle recovers P_i from the rest", and that this is not automatic for arbitrary two-term complexes. The elementary argument is therefore in the literature notes. What is new is the explicit loop hypothesis and the check on the walk data. The "Prior record" section does not cite this note and should. The caveat in E-161/E-166/E-167 and `research/syntheses/001-rounds-001-052.md:86` ("generation assumed") is the thing closed. The submission is honest that it is a caveat-closing, not a new phenomenon.

## Evidenced?
Yes for what is claimed. Ranges, counts and the failing-dimension breakdown are stated and reproduce. The failing-step set depends on a 519 s run, and I did not independently rebuild it.

## Scope
The title says "for every loopless step (so for all tiltingPlus steps on acyclic algebras)". Narrow it to "for every step whose mutated vertex has no loop", plus a statement that loop-freedom of children is checked only on the 29 edges. "Acyclic algebras" does not imply that the children stay acyclic.

## Required for acceptance
1. Cite the rickard-morita note as prior art, and rephrase "new" accordingly.
2. Either check that every child on the E-163, E-160 and E-157 paths has no loop at the next mutated vertex, or restate the corollary conditionally.
3. State that the 2116 steps are not gate-admitted, and say that H1-H3 are premises of the proof and are checked only for confirmation.
4. Note that the empty-M case is also covered.
All are doable in one sitting.
