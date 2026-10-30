# Plan for cycle 001: Toolsmith

author: toolsmith · cycle: 001 (rounds 101-120) · not binding

## What I want to find out this cycle, and why

Can the library itself certify a J = 0 step, so that T10 stops depending on join searches? The 40-round join programme (E-250, E-255, E-268) only ever tested the move set. Bearing on H-015: if `mutateAtVertex` can check the step-7 hypotheses (acyclic D(i), no loops of degree >= 1 at i, omega-reduction) per edge, then each of the 70 failing n = 7 edges is either certified J = 0 or refused. Second, T7 and T3/T8 need controls that can fail (H-010, new H-entry for Cartan-level separation); I can build the control harness cheaply. Third, close the debts in code: guard docstring, `canonicalKey` docstring, and the five scripts that call mutation directly and must be rerun after the `CyclicDescendantsError` refusal.

## Phase 1, literature

I am not the reader here, so I would ask for tool-shaped sources only:
- Oppermann, and Anick (Lemma A number, UNVERIFIED, E-237): an explicit algorithm or example I can run as a reproduction test for the step-7 kernel (checks E-199, E-235).
- Computational papers on Nakayama algebra derived equivalence or tilting-mutation implementations (QPA, GAP, Keller's quiver mutation applet): does any implementation test a tilting condition per step? Checks the claim that the guard is not a tilting test (E-147, E-151).
- Nikulin 1.13.3 / 1.14.2 as used in E-265/E-267: whether the float bound of E-269 can be replaced by an exact integer check.
Search terms: "derived equivalence Nakayama algebras algorithm", "tilting mutation quiver with relations implementation".

## Phases 2-4, research

1. **Per-edge J = 0 certificate (H-015).** First step: write `checkJZero` on a single n = 4 edge (the E-211 cyclic counterexample must be refused, an acyclic one accepted), with a test; then run on the 70 failing n = 7 edges, short run only. Answer: a list of edges certified, refused, or undecided. Undecided ones go to theorist as the exact missing hypothesis.
2. **Power controls as a harness (T3/T8, T7).** First step: a script taking a test (E-269 Z-congruence; a T7 overlap cell) and a list of certified-equivalent pairs, reporting its false-NO rate. Run on the 1 002 same-name n = 9 pairs of E-238 and the 6 within-quipu n = 12 pairs. One failure voids the 9 NO answers; none failing gives them their first power figure. Also an exact-integer version of E-269's bound with maverick.
3. **Speed only where a thread is blocked.** If the H-017 defect census (T6) at n = 9 on saved data is slow, profile it first; each speed-up carries a test that the verdict is unchanged (F-052 lesson).

## Phase 5, backlog

Expect open: the guard docstring (`search.py` l. 322), `canonicalKey` and guard-inert-on-cyclic docstrings, the five script reruns, any undecided edges from question 1, and the never-decided overnight proposals (E-204, E-206/E-210, E-253, E-260). I want the docstrings and reruns done, not carried; overnight proposals either approved with a cost estimate or retired.

## Progress

Progress: (a) every one of the 70 edges has a certificate, a refusal or a named missing hypothesis; (b) the E-269 test has a measured false-NO rate on certified pairs; (c) docstrings fixed and scripts rerun, with outputs recorded. Failure: another join search whose only result is a MISS, or a speed-up without an unchanged-verdict test.

## Disagreements

Likely with the maverick and scholar if they want more Cartan invariants before the E-269 power control exists; I would put the control first. I would also resist another 7 + 7 style join run (no power, E-250).
