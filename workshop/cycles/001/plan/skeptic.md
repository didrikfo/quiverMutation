# Plan for cycle 001: skeptic

author: skeptic · cycle: 001 (rounds 101-120) · not binding

## What I want to find out this cycle, and why

H-015 (guard sufficiency, OPEN): can a key-guarded walk leave a derived class through a J != 0 step? Every "25 of 25 joined" rests on the J = 0 premise, and join searches cannot test it (E-250, E-255). I want one non-move-set certificate with power, or a plain statement that none exists. Second, the Cartan-level NOs (E-265, E-269) rest on an unsourced necessity step and a float bound, and no power control was run on certified-equivalent pairs. One failure voids nine NOs. I would try to break both before more is built on them. New hypothesis if wanted: "derived equivalent LNAs are Phi-equivariantly Z-congruent"; it fails if one certified-equivalent pair is not.

## Phase 1, literature

- Source for Phi-equivariant Z-congruence as necessary for derived equivalence (checks E-267, the unsourced step): Nikulin discriminant forms (1.13.3, 1.14.2 are UNVERIFIED), Euler form/Coxeter transformation invariance under derived equivalence (Happel, Ladkani), Z-congruence of Cartan matrices (Ladkani, Hille-Ploog style).
- Independent inequivalence invariants for Nakayama or gentle/gentle-like algebras (Hochschild cohomology, Avella-Alaminos-Geiss invariant, Bessenrodt-Holm), for T10/T5 and E-259's eight pairs.
- Tilting mutation of algebras with relations beyond Oppermann/Anick (Lemma A number, E-237; loops of degree >= 1 at i, T5).
- Pull the LaTeX and match theorem numbers before citing; mark anything else UNVERIFIED.

## Phases 2-4, research

1. T3/T8 power control (E-269 test on the 1 002 n = 9 same-name pairs of E-238 and the 6 n = 12 within-quipu pairs). Pre-register: any non-congruent pair is a refutation of the necessity step. A clean pass only shows the test does not misfire on equivalent pairs, and I would say so.
2. T5/T10: build a certificate for the E-229 / 70 key-kept J != 0 edges that does not use the move set (silting engine, Hom(T,T[m]), E-087, E-053). Hit means a child outside the class. A miss on all 70 supports the premise, but only at n = 7.
3. If time: the E-258 open cases in `skeptic_convention2.py` (non-monomial, cyclic, parallel-arrow parents) and E-235 step 4 against Lemma A's hypotheses.

## Phase 5, backlog

Expect left open: the E-235 step-4 gap, S-1's K-threshold law (three lengths only), and E-151 counts that are capped BFS samples. I want the `mutationSearchDepthFirst` and `canonicalKey` docstrings corrected, UNVERIFIED citations resolved or struck, and the five toolsmith reruns checked.

## Progress

Progress: a Cartan-level NO with a real equivalent control, or a documented refutation of it; one F-entry for the J = 0 theory (cancellation, Lemma A, F/R convention); a certificate on the 70 edges. Failure: another join search, or another invariant that has no power on the group it was tried on.

## Disagreements

I expect the maverick to want breadth on the H-017 defect census. I would put it after the power control, since saved n = 9 data is cheap but proves nothing about the premise.
