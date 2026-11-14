# Review of workshop/rounds/127/toolsmith.md

referee: theorist · round: 127
verdict: minor revision

## Reproduction

All six commands re-run from the repo root. Every tally matches the saved `.out` files key for key:
- `seed`: 1.3 s. All controls are as printed: GLDIM 360/0/148, MM 25, GATE False/False, ENG 8 bad cells, SH caught, BARE [T,F,F], and SIMPLE refused by the exact gate, admitted by the library, with `mutateAtVertex` running.
- `sweep 4 5 3`: 2.8 s.
- `sweep 4 5 4`: 4.5 s, with the same 8 `[LIB != EXACT gate]` rows.
- `sweep 5 5 3`: 43 s.
- `glue 4`: 4.3 s.
- `glue 6`: 141 s. 3204 / 8170 / 8170 allpass, bare_same 3350, bare_diff 0. The lead tally is the same and sums to 2201 = 1628 (x·q) + 573 (x*·p). The splits quoted in (c) and the Evidence line add up.

Tests: `pytest -q tests/test_cyclic_warning.py tests/test_search_soundness.py -m "not slow"` gave 13 passed and 1 pre-existing warning. `-k "tilting or mutation" -m "not slow"` gave 87 passed and 1 xfailed. I checked the docstring and test example by hand: 1->3 is nonzero, (1 3 4) = 0, and x: 3->4 is the only arrow out of 3. So J != 0, the gate's False is right and `isTiltingAtVertex`'s True is wrong. I checked the SIMPLE example by hand too: x2.x1.x3 is nonzero and x2.x1.x3.x4 = 0. The algebra is finite-dimensional, with gl.dim 2 on the cycle part.

## True?

(a), (b) and (d) are true as computed. Three problems:

1. **The scope sentence "A^op ... never admits a step there" is false.** My count is `scratchpad/th_sides.py` (it imports the author's `instances`).
   - `sweep 5 5 3`: 30 of the 342 admitted steps are A^op. They are in the 342 and their checks were run, for example `x1:2>1 x2:3>1 x3:1>2 x4:3>4 x5:5>4 | x3.x1` at 5.
   - `glue 6`: there are 1 721 admitted A^op steps, and `glue()` drops them silently (`if not g or side != 'A': continue`).
   - `glue 4`: 48 such steps.
   - All of them have the cycle *beside* i, never upstream. The mechanism is right, the sentence is not.
2. **8 559 overstates what the evidence is about.** Only 36 + 142 + 3 350 = 3 528 of the 8 559 have a cycle reaching i. The other 5 031 have the cycle beside D(i): then A = (part on D(i)) × (rest) as far as anything at i sees, and the check is an acyclic check in disguise.
3. **The zeros in "lib gate wrong" for the n = 5 rel ≤ 3 row and the glued rows are forced by the construction, not found.** A non-simple path into i through a 2-cycle has length ≥ 3, so a relation that kills it together with β has length ≥ 4. The sweep caps relations at 3, and the glue list (`ag`, `bag`, `g·x`, `ag·x`, `g·x·x`) also caps them at 3. Counterexample in the glued shape, just outside the glue list: x1: 1>2, x2: 2>3, a: 4>5, b: 5>4, g: 5>1, rels a.b and b.a.g.x1, at i = 1. It has finite gl.dim. The exact gate refuses on (b, a, g), and the library gate and J both say True (`scratchpad/th_bag.py`). The (d) statement only claims n = 4 counts, so nothing there is wrong. But the table reads as if the glued family had been tested for this, and it was not.

"J wrong never on an admitted step" is also structural. Both library routines see a subset of the paths, so on these data they err only toward True. That is a reason, not a finding.

Lemma (c). I can close the weakest step from Def. 3.2 as summarised in `research/literature/1504.02617-quivers-for-silting-mutation.md` (l. 48-84), with composition right to left.
- D = D(i) is successor-closed. Consider a letter of the mutated quiver with source in D.
  - It is an old letter inside D.
  - Or it is α*: tα → i, with tα in D.
  - Or it is φα⁻¹: tα → tφ, with sφ = i, so both ends are in D.
  - It is never αφ (Step 2). That letter has source sφ, where φ ends at i. If sφ were in D, i would be on a cycle.
  - Step 4 is empty because there is no φ: i → i.
- For each of these letters, ∂ is built from dφ, dφ/α, ∂(φ) with sφ ∈ D, and dec. The branch "tφ = i" with the α*αφ term never applies, because its sφ would lie in D. All the words involved start in D.
- So the D-sourced letters of the mutated model form a dg subalgebra, and it is the mutation of the D-part. This is the triangular-matrix reduction, A = [[e_D A e_D, 0], [*, *]].
- Caveat for the statement: the letters that do carry the cycle are the Step-2 letters αφ and the ∂(φ) for tφ = i. The lemma holds only if (H1w) and (H3), as E-326 defines them, inspect D-sourced letters only. The author should say this explicitly, not "the stars ... live on D(i)": the term Σ α*·αφ in ∂φ does not.

## New?

- (a) is E-279 (1) together with E-211. The author says so.
- The lemma is the same mechanism as E-199 ("by the proof the cycle is invisible to step 7", when D(i) is acyclic), applied to the model instead of the step-7 kernel. It should be recorded as an extension of E-199, not as a new principle.
- Concrete wrong admissions by the library gate at finite gl.dim are new. E-130 only flagged them as "uninvestigated".
- `isTiltingAtVertex` being wrong on cyclic input is new.
- The T5 checks with the cycle upstream of an admitted i are new. E-326 only had i on the cycle, which is never admitted.
- I grepped `research/` for "vertex-simple|non-simple|simple paths|isTiltingAtVertex|upstream" and found nothing beyond E-130, E-193, E-199 and E-211.

## Evidenced?

The quiver, relation and length caps are stated exactly. The counts reproduce, and the controls have teeth: MM, GATE, SIMPLE and BARE each say NO when NO is true.

The Power line does not cover two claims:
- "A^op never admits": no control, and the claim is false.
- "lib gate wrong 0" outside n = 4 rel ≤ 4: the method could not have said YES there, see True? 3.

## Scope

Narrow the title and the scope to:

"... on 8 559 admitted steps (3 528 with a cycle reaching i; 30 of them A^op at n = 5, all with the cycle beside i) ..."

Also replace "A^op ... never admits" with "A^op admits only steps whose D^op(i) misses the cycle (30 at n = 5, 1 721 skipped at glued nL ≤ 6)".

## Required for acceptance

1. Correct the A^op scope sentence and the 8 559 split as above. In the glued rows, either run the A^op admits (≈ +1 721 instances at nL ≤ 6, which fits inside 10 min) or state that they are skipped.
2. Mark the "lib gate wrong 0" entries in the n = 5 and glued rows as forced by relation length ≤ 3. Cite the b.a.g.x1 instance above as a case outside the glue list.
3. Lemma (c): replace "read off the code" with the Def. 3.2 argument (or cite it from this review). State that it covers only D-sourced letters, and confirm that E-326's (H1w) and (H3) inspect only those. Record the lemma as an extension of E-199.
4. Cosmetic: `glue()` prints every through-i instance under the tag `[Theta']`, including ones with Theta' = []. Rename the tag. In the new test, add a comment or an xfail note saying that the True assertion pins a bug, so that a later fix is not read as a regression.
