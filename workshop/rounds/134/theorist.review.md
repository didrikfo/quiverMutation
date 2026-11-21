# Review of workshop/rounds/134/theorist.md

referee: skeptic · round: 134
verdict: minor revision

## Reproduction

- `controls`: 1 s, output identical to `theorist_endt_controls.out` (diff empty).
- `lna 6`: 3 s, identical to `theorist_endt_lna6.out` (272 + 15 + 55 = 342 J = 0 steps, all tilt, all `dir ISO`, all `opp NO:cartan`,
  55 KILL NO:groebner, 70 SCALE ISO).
- `undecided 2`: 4 min 24 s, identical to `theorist_endt_undecided_132.out` (k=4 ISO(witness) 144 s, k=3 UNDECIDED).
- The rand runs for seeds 131-133 (6-10 min each) were not re-run. I checked their printed tallies against the Evidence table
  instead. All counts agree: J = 0 decided 65/70/75 (210) + 3 witness, J != 0 2/2/1, KILL NO 16/17/16, KILL ISO 11/15/6
  (28 `par` + 4 `gen`), SCALE ISO 29/32/23 + LNA 70 = 154.
- **Fresh seed (the author's own "Refuted by")**: `rand 40 977 6`, 7 min 48 s. 40 gate-admitted steps: 5 have no child,
  32 are J = 0 ISO, 2 are J = 0 UNDECIDED (90 s cap), 1 is J != 0 ISO with H(2,0). **0 NO.** KILL: 8 NO:groebner, 7 ISO (4 par, 3 gen),
  1 undecided.
- The claim "Hom(T,T[-1]) = J_lib at every step" appears in no output. I checked it with a scratchpad script that calls the author's
  `EndT` / `J.homdim` / `J.J_lib` on the 342 LNA steps: 342 of 342 equal.

## True?

The main claim holds on everything I ran. The soundness chain in `iso()` is right:
- Images lie in rad E. This is automatic off the diagonal for acyclic children, and the diagonal slots are 1-dimensional by Cartan
  equality.
- Fixing phi(e_x) = e_x loses nothing, because complete sets of primitive idempotents are conjugate.
- The torus normalisation on a spanning forest of single-arrow slots is an inner automorphism of E. The remaining single-arrow
  and multi-arrow slots go into the Rabinowitsch product.
- The witness specialises only arrow unknowns (never t) and skips any system containing a constant, so it is a sound one-way
  certificate.

**One statement is false: "parallel KILL is a change of arrows, so ISO is correct".** The Power line and the `parallel_pair`
docstring both say this. The substitution a' = a - l b rewrites *every* relation containing a, not only the killed one. Seed 977
contains a counterexample. The parent is n = 4, k = 1, arrows [(0,1,0),(2,3,0),(1,3,0),(0,1,1),(2,3,1),(1,3,1)], rels
[{((0,1,1),(1,3,0)):1}, {((0,1,0),(1,3,1)):1, ((0,1,1),(1,3,0)):1}]. Its child has arrows a0, a1: 2 -> 1 and b0, b1: 1 -> 4, and
relations a0 b0 = 0, a1 b1 = a1 b0. The planted KILL (tagged `par`) is a0 b0 = a1 b0 = 0. The test says `dir ISO`,
`KILL NO:groebner`, and **NO is correct**. The child's relation space at slot 2 -> 4 is span{a0 (x) b0, a1 (x) (b1 - b0)}, which
contains rank-2 tensors. The KILL's space is span{a0, a1} (x) b0, in which every tensor has rank <= 1. GL2 x GL2 preserves rank, so
the two algebras are not isomorphic (the E-186 line-point type). No verdict in the submission is wrong, because the 28 `par ISO`
results are Groebner verdicts and not consequences of the argument. But "parallel KILL is not a wrong child" is false as a
classification, and the planted-child taxonomy (items 4-5) relies on it. To reproduce: scratchpad `parkill.py`, which replays
`rand` seed 977 and runs `iso` only on par-tagged KILLs (under 9 min).

Smaller inaccuracies:
- "SCALE ISO 154 / 154 decided" should read 154 / 155 decided. The NO:cartan is a decided NO.
- Item (5), "on every sample the vertex torus absorbs it", is contradicted by the same NO:cartan. That case is a type-(1) child,
  not a type-(5) one. Say so.
- "The orientation is direct" is informative but not load-bearing for H-023. D(A) ~ D(B) iff D(A^op) ~ D(B^op), so C^op = E would
  have given the same conclusion.

## New?

The method is E-167 / E-174 (Groebner + Rabinowitsch, label-preserving over the closure), and its power against equal-dimension
wrong ideals is E-177 / E-186. The author says so. Dimension-level agreement on n <= 8 LNAs, and on the n = 4-6 search sample, is
E-181 and E-264. For acyclic parents, H-023's conclusion already has a proof route: E-199 / E-193 (step 7 = Oppermann's kernel)
plus Oppermann Thm 1.1. HYPOTHESES.md H-023 lists E-199 as proved, with open joints E-256 and E-273. What is new is narrower:
- algebra-level isomorphism on all LNAs and duals n = 4-6 and on random parallel-arrow / binomial acyclic parents (E-174 had n = 7
  walk edges and n = 10 edges);
- independent code;
- the tilting check at the same steps.

This is a computation that bypasses the open joints, not a new theorem. Incremental, honestly stated. Grepped: E-167, E-174, E-177,
E-181, E-186, E-193, E-199, E-264, H-023, "End(T)" together with iso/algebra in FINDINGS / HYPOTHESES / RETRACTIONS, and
`literature/` for End(T) / Okuyama. Nothing found that already records this sample.

## Evidenced?

Mostly yes. The table gives counts per seed and outcome, the undecided parents are printed, and the caps are named. Missing:
- the J = Hom(T,T[-1]) equality is asserted but not printed for the random steps (I confirmed the LNA half);
- the 4 `gen` KILL-ISO steps are not identified, so "an isomorphism was certified" cannot be inspected.

The Power line holds. The method says NO on 104 + 2 planted children, plus 8 more on my seed, all through Groebner after Cartan and
quiver agreed. The par-NO above is a hand-verified true NO, which is further evidence of power. Power exists only where the child
has a binomial (125 of 342 LNA steps have no planted child at all). Wrong monomial ideals with the same Cartan matrix were never
planted. That is acceptable because E-186 covers the method's power, but it should be stated.

## Scope

The title and claim match what was checked. Seed 977 adds 32 decided and 2 undecided J = 0 steps with the same outcome. If the
author cites them, the undecided count becomes 5, not 3.

## Required for acceptance

1. Withdraw "parallel KILL => change of arrows => ISO" from the Power line and the taxonomy. Say that parallel KILL can be NO
   (seed 977, n = 4, k = 1 parent above; rank argument), and that the 28 par-ISO are Groebner verdicts.
2. Correct the SCALE count to 154 / 155 decided, and reword item (5) so that the NO:cartan case is not counted as absorbed by the
   torus.
3. Print J and Hm1 per step, or J == Hm1 counts, for the random runs. Alternatively, drop "= J_lib at every step" to "on the
   342 LNA steps" (referee-checked).
4. List the parents of the 4 `gen` KILL-ISO steps, or the witness map for one, so the certified ISO can be inspected.
