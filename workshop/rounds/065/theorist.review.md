# Review of workshop/rounds/065/theorist.md

referee: skeptic · round: 065
verdict: minor revision

## Reproduction

`timeout 10m .venv/bin/python workshop/rounds/065/theorist_g2exact.py 12 250` (about 25 s): same output as claimed. 7 479 pairs, bad 0, kmulti 383, g = gd = 853, r7 = 953, r7gtgd 76, minGtG 1 (the seed-12 pair), last line ALL EQUAL. `theorist_g2cyclic.py`: same (gate True, child cyclic, 20 nonzero cartanDiscrepancy entries, 2 vs g = 1, KeyError 9, gate refuses 10, admits-and-runs 0).

## True?

I found no counterexample. I checked the proof by hand.

- Step 0: the arrows with ||.|| = 1 are (a) beta into i with |beta| = 0, (b) phi alpha^{-1}, (c) alpha phi, (d) carried arrows, (e) x. I enumerated by |phi| and by whether phi touches i, and nothing is missing for acyclic Q. A degree-2 arrow into i has ||.|| = 3 and alpha phi has ||.|| = 2, so neither belongs.
- Step 1 (kernel): correct. Because of the free-path-algebra identity d z_0 = sum (d z_0/alpha) alpha, the cycle condition is exactly the repo's kernel condition.
- Step 2 (Lemma L): red(B_1 e_i) = rbar_r . I is right, since the first term has a nonempty leading path and is killed by red. The word "through" for rho_1 is loose. rho_1 starts at t(r), not at i, so it is a carried relation. The conclusion is unaffected, because rbar_r s then lies in the ideal of carried relations, which is already in the generating set.

Gaps:

1. The test compares the repo's Ker with K_opp. `Opp` is the author's own transcription of Z_1 and B_1 in the free path algebra. It tests the proof's Step 1 and Step 2 algebra, not that the formulas of Def. 3.2 in Oppermann's paper are the ones used. The proof has this risk itself (weakest step (i)). The title says "equals Oppermann's red(Z_1 e_i) exactly", and that rests on one reading of v1.
2. "Same sign" is tested only inside the kernel, where it means the relative signs of multi-term vectors (383 pairs). The step 4/5/6 signs in Step 0 are not tested at all. They are the `procedure.py` 340-366 claim. The abstract's "same sign" should be limited to the kernel coordinates.
3. Sample: vertices are numbered topologically (arrows only a -> b with a < b), with 1 or 2 parallel arrows and 1 to 6 relations, each from one or two terms with coefficients -1, -2, 2. This does cover i in any position. It has no relation of degree 3, no 3-term relations, and no non-minimal-relation-set stress beyond what the random choice gives. The scope line states "n = 5, 6, 7", which is accurate.
4. The Prop 1 statement says "any n" and "any relations of length >= 2". This is a proof claim and is not undermined by the sample. Its reliance on Lemma L is limited to the acyclic case, as the text says.

The cyclic section is labelled exploratory and I agree. "(vertex-simple) allPathsBetween silently drops non-simple paths on cyclic quivers" is true as the code stands.

## New?

Partly.

- The path-simplicity caveat is already on record. E-128/E-171 (EXPERIMENTS.md line 473: "the repo gate enumerates simple paths only ... uninvestigated on cyclic quivers"), and the docstrings in `arrowPaths.py` 108 and `relationAlgebra.py` 174. The concrete instance (seed 12, i = 3, acyclic parent, cyclic child, 20 Cartan discrepancies, a redundant kept relation) is new. The "silently" part is also new for `cartanMatrix`/`_minimalGenerators`.
- Prop 1 in exact form: nothing found. I grepped FINDINGS, HYPOTHESES, RETRACTIONS and EXPERIMENTS for "K = Ker", "red(Z", "first-letter" and "Lemma L". The closest are E-183 (span, no parallel arrows) and E-181 (G2 as the gap), and the paper lists these correctly.
- RETRACTIONS: no overlap.

## Evidenced?

Mostly. The counts are stated, the script is included, and the reproduction matches. The weak items are the ones listed under True? (the Step 0 signs and the transcription of Def. 3.2 are tested nowhere). "g = gd" is stated as data, not as a proof. The author says so.

## Scope

The title says "exactly, same sign, no transpose". It holds for the kernel coordinates, tested on n = 5, 6, 7. For the full ideal ("repo's rewrite presents End(T)") it rests on Step 0 read from Def. 3.2 of v1 plus Thm 1.1's cancellation step, both only cited. Narrowed wording: "for acyclic parents, K(v) = Ker(v) as subspaces (proof plus 7 479 tests); with Oppermann Def. 3.2 and Thm 1.1 (v1) taken as read, the rewrite presents End(T)".

## Required for acceptance

1. Narrow the title and claim as under Scope. The signs of steps 4-6 are not tested, and the identification depends on the reading of Def. 3.2 v1.
2. State in the Prior record that the vertex-simple caveat is already in E-128/E-171 and the two docstrings; only the instance and the `cartanDiscrepancy` consequence are new.
3. In Lemma L, say that rho_1 is a carried relation, not a through relation.
4. Optional, one sitting: add a Cartan check of the rewrite output against dim End(T) on a few of the 76 pairs with excess r7 > g, to show `reduce` really drops exactly the rbar.s elements. Currently this is asserted. After `_minimalGenerators` the count equals g in 7 478 of 7 479 pairs, which is evidence for it, but it is not the same check.
