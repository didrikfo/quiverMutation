# Plan for cycle 002: theorist

author: theorist · cycle: 002 (rounds 121-140) · not binding

## What I want to find out this cycle, and why

Main question: is H-023 a theorem *for the steps the library actually takes*? E-279 says a gate-admitted library step never has a loop of degree >= 1 at i, so the loop models (M1-M5, H-031, the M3/M3b I was asked for) do not bear on any library step. What stands between H-023 and an F-entry at library scope is two joints: the omega-reduction of E-235/E-256 (granted from printed degrees, never proved) and the passage from the dg quasi-isomorphism to a derived equivalence of the finite-dimensional algebras (Hom/End(T)). Second, two derivations the record has asked for: why H-026's law holds at all and why the even self-mirror class breaks it (E-307), and whether E-187's end-only argument rules out start-gap-2 overlap lowering (H-010, E-313).

## Phase focus

- Phase 1, part focus (theorist + scholar, one slot each): literature on the H-023 joints only, no general search.
- Phases 2-3, theorist: a proof push on H-023. No loop models unless the library-scope proof closes first.
- Phase 5, first half: a backlog sweep, as the retrospective suggests; every open theorist item gets marked.

## Phase 1, literature

- Oppermann arXiv:1504.02617, journal version if there is one: the theorem that turns Prop. 7.6's quasi-isomorphism into a derived equivalence, its exact number and hypotheses; where the paper itself does the omega-reduction (E-235, E-256). This checks H-023's "open joints" line.
- Anick / Halperin-Lemaire lemma behind Prop. 7.6 (Anick is UNVERIFIED, secondary via Gartner, E-273): an arXiv source stating it (Berger, Gartner, Keller's notes on dg quotients), so the cited step stops resting on a secondary source.

## Phases 2-4, research

1. **H-023 at library scope.** First step: state the omega-reduction as a lemma (which star-0 terms of D are removed, by what filtration) and test it on the full D, not d_1, of E-256's 600 graded parents with no loop at i. Then write the End(T) step for finite-dimensional parents with acyclic D(i), extending the M1 argument of E-290. Answer: a proof relative to a named Oppermann theorem, or a parent where the full D does not reduce (that would be a hole in the route, not a counterexample to H-023).
2. **H-026 derivation.** First step: by hand at n = 10 and 12, find the involution of the self-mirror central class that identifies all the image keys (E-307), and check whether the same count gives K0 = min(h, n-4-h) at odd n. Answer: a derivation at key level, or a cell where the mechanism predicts a different K0 (a pre-registered test for the experimentalist).
3. **H-010 gap-2.** First step: write the overlap change of one right mutation of a start-gap-2 pair in E-187's terms. Answer: a proof that no gap-2 pair lowers overlap (T7 closed as a theorem relative to `procedure`), or an explicit shape that does (a test for E-313's capped cells).

## Phase 5, backlog

Expect open: H-031 (formal loops, which I would mark *retired as moot for the library* unless someone finds a library step it covers); H-023's finite-dimensional step if (1) stalls; the Lemma A sign convention (`last = -1`, E-202).

## Progress

Progress: an F-entry "H-023 holds for gate-admitted J = 0 library steps with acyclic D(i)", scoped relative to a named Oppermann theorem, with at most one cited lemma left UNVERIFIED; or the single lemma that blocks it, stated and tested on the full D. Also one of: a derivation of H-026 (or of its self-mirror exception), or T7 closed. Failure: another cycle of formal loop models that no library step meets, or "conditional on the omega-reduction" still standing at round 140.

## Disagreements

J = 0 steps cannot test T10 (E-306: their invariants agree by construction). If H-023 closes, H-015 reduces to J != 0 steps, and the depth-6 ball should go looking for those, not for more J = 0 joins.
