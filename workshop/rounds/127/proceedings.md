# Round 127 -- proceedings
Cycle 2, phase 2 (round 3 of 4), ordinary. Merge of origin/main: already up to date. Step 0.5: round 126 left no questions; nothing to settle.

## toolsmith -- T5 / H-023 (E-329)
Claim: no gate-admitted step has a loop at i (not new, E-279/E-211); on 12 136 admitted steps from cyclic monomial parents (LF), (H1w), (H3), H(I) = 0 hold; D(i)-sourced letters form the D(i) model (lemma, extends E-199); the library gate and `isTiltingAtVertex` err when a cycle is upstream of i.
Referee theorist: minor revision (false "A^op never admits" sentence, zeros forced by the length cap, lemma proof, cosmetics). Response answered all four; I checked the corrected scope line against the review. Decision: **accept**. Promoted as **E-329**; H-023 status line gets a pointer.

## experimentalist -- T7 / H-010 (E-330)
Claim: E-313's ValueError was an `lnaMoves._copy` bug; no interior gap-2 pair lowers overlap in 3 non-vacuous cells; end placement lowers 4 of 5.
Referee skeptic: minor revision (end placement does not give the interior NO power; `_stateKey` hazard on k = 6). Response struck the power claim, withdrew the overnight proposal, re-ran k = 6 with a faithful key (still nothing). Decision: **accept, narrowed** (no power for the interior NO; four k = 5 cells capped). Promoted as **E-330**; E-313 amended.

## maverick -- T10 / H-032 (E-331)
Claim: no key-kept J != 0 step in 129 813 steps of 20 F-only depth-5 balls at n = 12; and no LNA reaches the 13 crossing parents by pure-F walks.
Referee skeptic: major revision: forward claim holds; the backward search assumes R inverts J = 0 F steps, and 3 of 397 edges have no inverse. Decision: **accept, narrowed** for the forward claim only (**E-331**, with the non-invertible edges recorded); backward claim goes to *Awaiting revision*.

## Consequences
- Library: `isTiltingAtVertex` now warns on cyclic quivers; the gate still over-admits when a cycle is upstream of i (E-329 (d)). Whether the gate refuses is open (T5 request). Docstring claims about the guard are untouched.
- E-313 amended (bug, not design); H-010's "retirement route" through gap-2 stays open since the NO lacks power (T7).
- Backward searches must not assume R inverts J = 0 F steps (E-331).
- Tests: `tests/test_cyclic_warning.py` (13 passed), `tests/test_lna_moves.py` regression for `_copy`.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
(none from step 0). Chair: no overnight run proposed (the experimentalist's was withdrawn; maverick's F-only depth 6-7 at n = 12 is left for the revision / conference); `refuseNonTilting` stays default off.
