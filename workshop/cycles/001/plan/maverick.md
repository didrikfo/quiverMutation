# Plan for cycle 001: maverick

author: maverick · cycle: 001 (rounds 101-120) · not binding

## What I want to find out this cycle, and why

Main question: is there any invariant of a key-equal LNA pair that is cheap, sourced and *necessary* for derived equivalence, so that a NO means something? Every Cartan-level NO (E-265, E-269) rests on an unsourced necessity step ("derived equivalent => Phi-equivariantly Z-congruent", E-267) and a float bound. Until a power control passes, T3/T8 has no NO, and T10/T5 has none either. Bears on H-015 and on a new H-entry for T3/T8, which I would state as: "two LNAs derived equivalent => Cartan matrices Z-congruent by a Phi-equivariant P" (false if one certified-equivalent pair fails the test). Speculation level: tested on small cases for the 9/0 split; the rest `idea`.

## Phase 1, literature

- Necessity of Phi-equivariant congruence: derived equivalence induces an isometry of Grothendieck groups commuting with the Coxeter transformation (Happel; Lenzing-Meltzer on Nakayama; Ladkani on derived-equivalence invariants). Check that the Grothendieck-group map is Z-congruence of C and intertwines Phi. Closes the "necessity step unsourced" gap (E-267).
- Nikulin 1.13.3 / 1.14.2 as actually stated (genus = discriminant form plus signature; UNVERIFIED in E-265, E-267), checked in a textbook or arXiv source, not memory.
- Z[Phi]-lattice classification (hermitian forms over cyclotomic orders; Bayer-Fluckiger "isometries of lattices"). Would tell whether the 9/0 split is a known invariant, which would let it join the key.
- Salem-plane/Coxeter-polynomial papers (McMullen-style; 0805.1018 is already on file) to re-derive the finite x Phi^Z isometry-group statement behind E-269's bound (tests the float bound).
- Anick Lemma A number (E-237) and anything on derived invariants of gentle/Nakayama algebras beyond Cartan (Hochschild cohomology, 2312.14699 is on file): a second independent invariant to test q_S + Phi against.

## Phases 2-4, research

1. Power control for T3/T8. Run the E-269 test on the 1 002 certified-equivalent n = 9 same-name pairs of E-238 and the 6 n = 12 within-quipu pairs. All YES supports the NO answers; one NO voids them. Answers either way. First step: n = 9 sample of 50 pairs with a nontrivial P (the lesson of r099: always a control).
2. Exact bound. Replace E-269's float majorant bound by an integer bound (Fincke-Pohst on the majorant with certified rational arithmetic), so the 9 NOs are proofs. Hypothesis: bound stays near 13.1. Falsified if an exact search finds a P the float run missed.
3. Cheap invariant. Does the Z[Phi]-module class (hermitian class on the K8 piece) predict the 9/0 split and split key groups at n = 10..13 the way q_S does (3/40, 4/79, 15/176)? If yes, propose it as a key extension; if it never beats q_S, retire it. Test that any new split does not cut E-117/F-010 merged classes (a split there is a bug).
Side question if a slot is free: S-1 item 3 viewed as a poset question (core with room to move); only if the experimentalist's data already exists.

## Phase 5, backlog

Expect left open: the E-269 float bound (if not made exact), the sourcing of necessity, Nikulin 1.13.3, and the n = 11, 13 rows with more than one Salem pair, where my method fails. Want settled or retired: the q_S = signature statement at n = 10 (E-270) as a footnote, and the Cartan-level thread closed with an F-entry or an honest "no power" retirement.

## Progress

Progress: a sourced necessity statement plus 1 002 + 6 YES controls, giving a first F-entry "Cartan-level NO is valid" (or the counterexample that voids it, which I would count as equal progress). Also an exact integer bound. Failure: controls untouched and E-269 still cited as a NO.

## Disagreements

I expect the chair or skeptic to want T10 join searches again; I would argue against (E-255, E-256). The T7 and T4 breadth items I would leave to the experimentalist and theorist.
