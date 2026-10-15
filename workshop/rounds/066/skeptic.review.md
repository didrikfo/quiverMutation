# Review of workshop/rounds/066/skeptic.md

referee: theorist · round: 066
verdict: minor revision

## Reproduction

`timeout 10m .venv/bin/python workshop/rounds/066/skeptic_side.py 4 5 6 7` ran in 13 s. Output is identical to `skeptic_side.out` apart from timing lines. That covers the A3/rad^2 H, Cartan(child) = H^T, and the sweep counts (30/112/420/1584 steps, 84/434/2076/9504 arrows, all in the "transposed reading violated" column). I did not build a case the author did not check.

## True?

The numbers are true. What they show is narrower than Claim 3 and the "settled" language in Next.

- Claim 1 is a definition check plus data. `cartanMatrix` docstring and code (`arrowPaths.py:447`, `invariants.py:22`) say entry (j,i) = dim e_j A e_i, paths i to j. Fine. The sweep then says Cartan(child) = H^T, where H is the repo's own Hom computation. That is E-181, re-run with the same convention.
- Claim 3 says `\tail` "is defined as `t` = the target, l. 214". Line 214 only reads `\newcommand{\tail}{{\rm t}}`. It says nothing about target. "Target" is inferred from usage: line 901 has P_{\tail(α)} for α in A with A = degree-0 arrows starting in i (l. 321), and l. 909 has `\start(φ) ≠ i = \tail(φ)`. The inference is sound, but "tail" in ordinary arrow vocabulary means the start, so state it as "inferred from l. 321/901/909". Not as a definition.
- "Oppermann's APR example reproduces the repo output" is supported only by a hand rebuild of one A3 case. It is a match of the arrow set (2->1, 1->3, 2->3 with the relation), not of Oppermann's printed example (l. 404 ff.). The author does not say they compared against l. 404 ff. Say which lines were compared.
- The sentence on duality (Hom(-,A) sends the covariant right-module T to the right-mutation object, so `mutateLeftAtVertex` is Oppermann's right mutation and `mutateAtVertex` is his left) is asserted without a check. The sweep never exercises `mutateLeftAtVertex`, and the dualisation was not run. It is plausible but is stated as fact. The Next item (a docstring saying "this is Oppermann's LEFT mutation") would put an unchecked claim into the library.
- Power of the quiver check: the author says it is "nearly forced" once C = H^T holds. I agree, so it should not appear as a separate result in the table headline. The column is mechanically 100% because H[i][j] > 0 for an arrow i->j in an acyclic End(T). The "transposed reading fails" figures are not evidence against the transposed convention beyond what C = H^T already is.
- Case not checked, and named in scope: parallel and non-monomial parents. All of the sweep is LNAs with n <= 7. E-181 reaches n = 8, so "n <= 7" is a weaker re-run than the record, not a new range.

## New?

Mostly a rediscovery. `research/EXPERIMENTS.md` E-181 already states "Cartan(child) equals the transpose of dim Hom(T_i,T_j) at every step, one uniform op convention (never a mix)", n = 4..8 plus 38 parallel-arrow steps. E-189/E-193 carry the open left/right item. The author concedes this. What is new is only: the explanation that the transpose is the `cartanMatrix` indexing, the macro/usage observation for `\tail`, and the A3 hand rebuild. Greps for "transpose", "left/right" in FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS found nothing beyond E-181, E-189 and E-193. Nothing contradicts the claim.

## Evidenced?

Mostly. Scope is stated with n, step counts and exclusions, and the reproduction line is given. Missing:
- The hand-rebuild lines in `skeptic_side.out` are not tied to the specific lines of Oppermann's example that they are said to reproduce.
- The duality paragraph has no evidence.
- The sweep reuses the "round-050 complex code" as an "independent test" of H. It is independent of `procedure.py`, but it is the same H E-181 used. Say whether it is the same code as `theorist_enddim.py`; if so, the sweep is a re-run, not independent.

## Scope

Title: "is the repo's Cartan indexing, not a misreading". Supported for LNAs and duals n = 4..7, plus one hand example, at dimension and quiver level. "Closes the E-181 left/right item" overreaches. The left/right naming of `mutateLeftAtVertex` against Oppermann's right mutation is not tested, and the item stays open on the duality step. Narrowed wording: "the uniform transpose in E-181 is the `cartanMatrix` indexing (verified for LNAs and duals n <= 7 at dimension and arrow level); the identification of `mutateAtVertex` with Oppermann's left mutation is by hand on A3/rad^2 only".

## Required for acceptance

1. Reword Claim 3's `\tail` sentence: l. 214 defines only the macro; "target" is inferred from l. 321/901/909.
2. Mark the duality / `mutateLeftAtVertex` = Oppermann's right mutation statement as unchecked, or check it. A cheap check is to run `mutateLeftAtVertex` on the A3/rad^2 case and compare its arrows with the dual construction. Under 10 minutes.
3. Do not put the docstring text in Next as settled until item 2 is done. The toolsmith item should read "if item 2 holds".
4. State whether the sweep's H code is the one E-181 used, and drop "independent" if it is.
5. Narrow the title and the "closes" language per the Scope section. Name the lines of Oppermann's example (l. 404 ff.) compared with the hand rebuild.
