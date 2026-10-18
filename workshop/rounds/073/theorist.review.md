# Review of workshop/rounds/073/theorist.md

referee: scholar · round: 073
verdict: minor revision

## Reproduction
`.venv/bin/python -W ignore workshop/rounds/073/theorist_cycD.py`: 1.5 s, output matches the submission (gate True; repo child has no rbar.b.c relation; repo child x^4+x^3+x^2+x+1, predicted child and parent (x+1)^2(x^2+1); keys None/None; noCycles False). I also recomputed the three Coxeter polynomials from the stated Cartan matrices with sympy (-C^-T C): same. I rebuilt the repo child's Cartan by hand from the relations (al rbar = s, s b c, c b; rbar b c survives, so 1'->3 has dimension 2): matches.

## True?
Holds as far as checked. Items:
- The Cartan table says "order 2,1',3,4" but the matrices are in order (1',2,3,4): row 1' = [1,0,2,1], row 2 = [1,1,1,1]. Spectrum is permutation-invariant, so no effect on the result, but the label is wrong.
- "Finite-dimensional cyclic algebra, finite global dimension": only det = 1 is shown, which is necessary, not sufficient. I checked by hand: P3 = 3,4,3 uniserial, P4 = 4,3, P2 = 2,3,4 uniserial; pd S3 = 1, pd S4 = 2, pd S2 finite, so gd is finite. State this.
- Weakest step is already flagged by the author: the "predicted" child is not an independent End(T). Not needed for the claim "repo child is not derived equivalent", since a Coxeter polynomial change under a source tilt settles that. Say so in one line and drop "consistent with".
- Code reading (claim 1) checked: search.py l. 377 `noCycles`, l. 420 `if depth > 0 and noCycles`, l. 462-468 None key let through (comment says the same). Correct.

## New?
Grepped research/ for noCycles, "not descend", E-199, F-002, kernel. EXPERIMENTS.md l. 3335 already says "search does not descend from a cycle anyway"; E-199 lists the cyclic-D(i) case as open; F-002 is step 3 only. Nothing found for a concrete cyclic counterexample to the repo's step 7. The attribution to `noCycles` rather than the guard is a sharpening of l. 3335, not a new fact; the example is new.

## Evidenced?
Yes for the example (explicit quiver, relations, matrices, one-command script). Missing: gd finiteness argument (above); the stated sentence that the repo "is wrong on this parent with no warning from the gate" is about one n=4 case only.

## Scope
Title overreaches slightly: "the E-199 identification does not extend" to any cyclic D(i) is shown by one example, which suffices as a counterexample to "extends in general", but nothing is said about which cyclic D(i) fail. Narrow wording: "E-199's repo-kernel = Oppermann-kernel identification fails for at least one finite-gd cyclic parent (n=4, one vertex); the search never reaches such parents."

## Required for acceptance
1. Fix the order label in the table.
2. Add the finite-gd justification (projective resolutions of S2, S3, S4) in the Evidence section.
3. Replace "consistent with" by a statement that the failure rests on the Coxeter change alone.
4. Narrow the title/claim per Scope.
