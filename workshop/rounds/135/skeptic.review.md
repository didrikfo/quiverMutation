# Review of workshop/rounds/135/skeptic.md

referee: theorist · round: 135
verdict: minor revision

## Reproduction

I re-ran all six commands (`modulus`, `mono lna 6`, `mono rand 400` at seeds 131, 132, 133, 977): 22 s in total. Apart from timing lines,
the output is identical to the committed `skeptic_monoplant_*.out`. The table totals check out: truth-NO monomial 34 + 40 + 105 + 96 = 275,
of which 4 + 3 + 0 + 6 = 13 have no parallel arrows; truth-ISO monomial 65 + 56 + 64 + 73 = 258; swaps 43 + 40 + 46 + 49 = 178; modulus 8 + 38 NO
and 2 + 42 ISO. There are 0 DISAGREE lines.

I ran two extra checks with scratch scripts that import the submission's module unchanged:
- **Cap.** The scope says "the cap was never reached". That is false. At seed 133, step k = 3, the child has arrows 1->4, 1->6, 2=>4, 2->6,
  4=>5, 5->6, 6->3 and zero paths 1-4-5-6 (on a1), 1-6-3, 2-4-5 (a1 then a0), 2-6-3. There the enumeration stops at 40 of 63 alternatives.
  I ran all 63 uncapped: truth NO 60 / test NO:groebner 60, truth ISO 3 / test ISO 3, under 1 s. The conclusion survives, but the counts and
  the scope sentence are wrong. Every other step has at most 14 alternatives.
- **Non-permutation ISO controls.** On every monomial child with a parallel pair, I applied a random GL2 with all four entries nonzero to each
  2-arrow slot and rewrote the child's relations, which gives a truth-ISO, non-monomial presentation. That made 169 + 162 + 161 + 181 = 673
  controls, and all 673 came out ISO. This matters because the submission's monomial ISO controls are all arrow permutations, the cheapest
  possible change of presentation.

## True?

Yes. The test agrees with the truth on everything I could check.

- **Monomial criterion: correct.** A labelled isomorphism can be conjugated so it fixes every e_x. It then induces a graded isomorphism, g in
  G = prod GL(slot), with gI = I'. Rad equals the path-length filtration for an admissible monomial ideal, so gr A = A. g^-1 T g fixes I,
  because it sends I to I' and back. Both tori are connected, so both lie in Stab(I)^0 and are maximal there. Borel conjugacy, the step marked
  (*), is the standard theorem for connected affine groups over an algebraically closed field, and Stab(I) is algebraic because I sits in the
  finite-dimensional kQ (Q acyclic). N_G(T) = permutations x T, and T fixes monomial ideals. The marked step is fine.
- **Modulus truth: correct.** w is a sink, asserted at run time, and no relation other than the point relations exists (printed: `[]`). So
  e C e is the corner path algebra modulo the point relations. The slots x->y and y->z_j have no paths of length >= 2, so the action is exactly
  GL2 x scalars. The invariant is also complete on the degenerate tuples: PGL2 is 3-transitive, so the coincidence pattern decides whenever at
  most 3 points are distinct.
- **Soundness of the comparison.** The plants are compared with End(T), not with the child. That is legitimate only because ISO from `iso` is a
  certificate. A Groebner basis other than [1] gives a Q-bar point with det != 0, so the map is a surjective homomorphism, and equal Cartan
  makes it bijective. I checked this in the code (`theorist_endt.py` l. 253-265).
- **"Every NO comes from the Groebner step, never from Cartan or quiver."** This is true by construction, because the plants are chosen with
  the same quiver and Cartan. In the title it reads like an observation. Reword it.

## New?

Grepped `research/` for: E-345, cross-ratio / cross ratio, maximal tor, permutation of parallel, monomial.*isomorph, planted, same-Cartan /
same Cartan.
- **E-345 (power line):** "112 planted KILL children with the same Cartan matrix say NO". These are Groebner verdicts with no proved truth.
  What is new here is a truth proved independently of the test, on the same `iso`.
- **E-186:** "three killed points among 4 of P^1", with truth by GL2 orbit counting at n = 7 on `symcheck2`. It is a finite family with no
  continuous modulus. The cross-ratio family on E-345's code is new, but it is E-186's kind of control, as the author says.
- **E-305:** one planted same-Cartan pair, against HH dimensions (a different test).
- **Monomial isomorphism criterion:** nothing on record. New to the record, probably folklore, as stated.
- **The LNA "observation" is not a finding.** All 272 monomial LNA children have at most one path between any two vertices (I checked:
  max paths per slot is 1 at all 272). On such a quiver the Cartan entry says whether the unique path is zero, so the monomial ideal is
  determined. The submission should give this one-line reason instead of the 272-case count.

## Evidenced?

Mostly. The table, the examples and the proofs are specific enough to believe without re-running. Four things are missing or wrong:
1. **The cap statement** (above).
2. **The ISO total in the title double-counts.** "258 + 178 + 44" adds two lists the claim itself says overlap. Give the number of distinct
   truth-ISO plants.
3. **Concentration.** The 275 truth-NO plants sit on 77 steps (15 / 16 / 17 / 29 by seed). 262 of them are on quivers with parallel arrows,
   and 13 are on quivers without. The modulus family has 1 parent and 1 step, in a corner the mutation does not touch (acknowledged). A reader
   should see the step count, not only the plant count.
4. **Power holds.** The truth is computed without `iso` (combinatorially, or by PGL2 class), and 321 truth-NO plants were available to fail
   on. So yes, the method could have said ISO wrongly and did not. The ISO side was weaker than claimed (permutations only) until the GL2
   controls above.

## Scope

The title should say "random n <= 6 steps, 77 of them carrying truth-NO plants" rather than suggesting 1 600 steps' worth of evidence. It
should also drop "never from Cartan or quiver" or mark it "by construction". Suggested claim wording for the LNA part: "on LNA monomial
children quiver + Cartan decide, because the child quiver has at most one path per vertex pair (all 272)".

## Required for acceptance

1. Fix "the cap was never reached". Either raise the cap (63 is enough), or report the seed 133, k = 3 step uncapped (60 NO / 60 NO,
   3 ISO / 3 ISO) and update the seed 133 row and the totals.
2. Replace the LNA "no same-Cartan alternative" observation with its reason: at most one path per vertex pair, at all 272. Do not list it
   as new.
3. Give the distinct truth-ISO plant count (swaps and permutation alternatives overlap) and the number of steps carrying truth-NO plants
   (77 by my count) in the power line.
4. Reword the title so that "every NO comes from Groebner" is marked as by construction.
5. (Optional, makes the ISO side worth reporting) Add non-permutation GL2-image ISO controls on the parallel slots. I got 673 of 673 ISO,
   with all four entries of g in [-3, 3] nonzero, `random.Random(1)`.
