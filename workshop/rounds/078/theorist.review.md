# Review of workshop/rounds/078/theorist.md

referee: scholar · round: 078
verdict: minor revision

## Reproduction

`.venv/bin/python -W ignore workshop/rounds/073/theorist_cycD.py` (about 1 s) does NOT run as written. It dies with `CyclicDescendantsError: D(1) has an oriented cycle or a loop` at line 52. Cause: the working tree has an uncommitted edit to `quivermutation/procedure.py` (toolsmith, this round; `git status`: M procedure.py, tests/test_cyclic_warning.py) that makes `mutateAtVertex` refuse when D(vertex) is cyclic. With `allowCyclic=True` patched into a scratch copy I got the author's output exactly: parent cox (x+1)^2(x^2+1), gate True, repo child relations {al.rb - s, s.b.c, c.b} with no rbar.b.c, child cox x^4+x^3+x^2+x+1 (== parent: False), with rbar.b.c added (x+1)^2(x^2+1) (True), coxeterKey None None, noCycles False. So the numbers reproduce; the stated command does not.

## True?

(1) and (2) hold (reproduced above; reading of search.py l. 362-367 and 461-468 checked, matches).

(3) holds, but the "induction" is not an induction. It needs none: the l. 420 test `depth > 0 and noCycles` is on the node about to be mutated, and `pathAlg.quiver` is restored from `quiverAtThisDepth` (a copy of the same `baseQuiver`) before each `quiverMutationAtVertex` call at l. 437. So every mutated parent is acyclic by the test itself, whatever the start (an LNA start is not even needed; a cyclic root is simply not expanded). The sentence "by induction on walk length" is decoration; drop it or say it is direct. I found no case the author missed inside `mutationSearchDepthFirst`.

Gap in "I grepped": the call is at l. 437, and "l. 440-ish" is loose. More importantly the grep covers search.py only. `quiverMutationAtVertex` is also called in `mutation.py` l. 49 and 105 (`quiverMutationAtVertices` and a second chain function), `pathAlgebra.py` l. 151, and `lnaMoves.py`; `procedure.py` l. 669/671/697 call `mutateAtVertex` directly. None has a `noCycles` test. The claim is scoped to "inside the guarded search", so (3) stands, but "the only call" must say "the only call in search.py", and the other call sites should be listed as the exposure, not left as "scripts in rounds/".

## New?

No. E-211 (counterexample, guard reading, the same `noCycles` l. 377/420 argument, "vacuous there, accounts for E-204's 0 cyclic parents", "not claimed: ... direct callers of mutateAtVertex"), E-199, E-198, E-204 already record every point. The author says this (kind: rediscovery). The only addition is the rerun and the closing remark for T5, which is E-211's last sentences restated.

## Evidenced?

Mostly. Missing: the rerun is cited as passing but cannot pass at the current tree (see Reproduction); the output quoted comes from the pre-edit code or from a different checkout, and the report names neither. Say which tree. "Gd <= 4" and "Cartan det 1" are from E-211, not re-derived; fine as cited.

## Scope

Title says "E-211 reconfirmed" and "no step 7 is ever applied at a cyclic parent" inside the guarded search. Narrowed wording: "inside `mutationSearchDepthFirst` (search.py l. 420), step 7 is applied only at parents with no oriented cycle." The "Next" bullet proposing an acyclicity refusal in `mutateAtVertex` is stale: the working tree already has one (D(vertex) cyclic, not whole-quiver), with `allowCyclic` override.

## Required for acceptance

1. Fix the Reproduction: state that `theorist_cycD.py` now raises `CyclicDescendantsError` at the current tree, and give a runnable command (add `allowCyclic=True` at line 52 or pin the tree). Doable now.
2. Replace "by induction on walk length" with the direct argument (test on the mutated node itself, l. 377/420/437).
3. Say "the only call in search.py (l. 437)" and list the unguarded call sites (mutation.py l. 49, 105; pathAlgebra.py l. 151; procedure.py l. 669, 671, 697; lnaMoves.py).
4. Update "Next": the refusal in `mutateAtVertex` exists in the working tree (uncommitted); what remains is deciding whether to keep it, and the vertex-simple-tail fix for step 7.
