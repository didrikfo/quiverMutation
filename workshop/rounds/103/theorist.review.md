# Review of workshop/rounds/103/theorist.md

referee: toolsmith · round: 103
verdict: minor revision

## Reproduction

Re-ran `theorist_loopdeg.py` (1.5 s) and `theorist_cycsweep.py 1 300` and `... 2 300` (seconds each). Output is identical to the figures in the submission. The E-211 parent gives the same gate, error and `KeyError` lines. The sweeps give the same counts: 495/394/101 and 491/401/90, and 160/145/15 and 181/161/20. I did not re-run the Keller `curl`; I fetched the same e-print separately.

## True?

I found no error. Nothing I tried contradicts the claims.

- **Theorem statements.** Oppermann Thm 1.1 (repo copy, l. 263-265) says "no loops of degree 0 attached to i". Step 4 of Def. 3.2 handles the cyclic case. Aihara-Iyama has "Any mutation of a silting subcategory is again a silting subcategory", and Thm "when silting is tilting" gives the criterion for tilting objects. Neither has a quiver hypothesis. Keller's Rickard theorem (`thm:Rickard`, l. 1047) needs only homology in degree 0. All three citations match.
- **Loop argument.** The induction is correct for Construction 2.2 models:
  - m = 1: d(phi) is a relation, which is a closed path at i of length >= 2.
  - m >= 2: d(phi) is a combination of paths made of lower-degree arrows, and each of those has a degree-0 path by induction.
  - So i lies on an oriented cycle of Q0, and then D(i) is cyclic and refused.

Weaknesses:

1. **The library check does not test claim (1).** "Default call refuses all 495/491 on-cycle cases" follows from the `CyclicDescendantsError` gate by construction. It would hold even if the loop lemma were false. The sweep supports only the counts and the `KeyError` rate. The lemma itself is not machine-checked: nobody built a Construction 2.2 model, listed its loops of degree >= 1, and confirmed that i is on a cycle each time.
2. **Parent class.** The sweep is monomial, n = 3 and 4 only. The task asked about parallel arrows and non-monomial parents. The proof does not need monomiality, but none of this was run, and the scope line does not say so.
3. **Iterated steps.** E-249 records that after one dg step the child can have degree-1 loops at other vertices, and the repo passes to H_0 at each step. The claim "the search never meets a loop at i" is true of the algebra the library holds. It says nothing about whether the H_0 child is the dg child, so it does not carry over to derived equivalence along a chain. This is outside the stated claim, but the title ("what is open is Lemma A for loops") reads as broader.
4. **Gate caveat.** The author's own remark that J is computed with vertex-simple paths on cyclic quivers (E-193) means "gate-admitted J = 0 with i on a cycle" is not a reliable class. That makes the falsifier in (3) vacuous for the library, as stated.

## New?

Mostly recorded. The following are already on file:

- E-273: Thm 1.1 excludes only degree-0 loops, Prop. 7.6, Lemma A is ours.
- E-235 and E-240: Lemma A, loops left open.
- E-256 a: d_1 shape of Lemma A.
- E-277: `allowCyclic=True` raises `KeyError`.
- E-249 and E-248: loops after a step.
- H-023.

Not found in `research/`: the argument that a loop of degree >= 1 at i forces i onto an oriented cycle (grepped "loop of degree", "degree >= 1", "loops of degree"; only E-248/E-249 turned up). That argument, the explicit failing term of Lemma A, and the sweep counts are new. The Keller/AI hypothesis matching adds little beyond E-273.

## Evidenced?

The theorem matches are specific and I confirmed them. The loop-to-cycle lemma is a proof sketch whose weak step ("d phi nonzero") is flagged by the author. The library claim is partly tautological (point 1 above). `Power: n/a` is acceptable for a literature item, but claim (1) is a universal statement about the library and rests on a hand proof only.

## Scope

Narrowed wording: "For a model of A built by Construction 2.2 (A finite-dimensional, basic), a loop of degree >= 1 at i implies i on an oriented cycle of Q0. Hence the library (default path and search) never reaches Step 4 of Def. 3.2 for the algebra it holds. Chains of steps through H_0 (E-249) are not covered."

## Required for acceptance

1. Add a machine check of the loop lemma. Build the minimal model for about 100 random cyclic and acyclic parents, including parallel arrows and commutativity relations, up to degree 2 or 3. Confirm that every loop of degree >= 1 has its vertex on a Q0 cycle. A check of the Step 3/4 quiver rule on a graded quiver (as in E-248) would do. [one sitting]
2. Say in the scope line that the sweep is monomial only, and that the default refusals are by construction rather than evidence for the loop lemma. [trivial]
3. State in the claim that the H_0 passage (E-249) is outside it, or say why it does not matter. [trivial]
4. Reword the title and claim as in Scope. [trivial]
