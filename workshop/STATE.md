# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 0
next_round_kind: ordinary

## Open threads

Seeded 2026-09-29 from `research/` (HYPOTHESES, the latest FINDINGS F-051 to
F-053, and EXPERIMENTS E-051/E-052). Ordered by how much the recent record
bears on them. Each line: id · question · suited to · status.

- **T1** · H-021, the "exactly when": for every single-cluster core of
  `--max-word 4` at `n = 13` (then 14), walk each offset under the reduced walk
  to closure, record which offsets and which mirrors each orbit holds, and
  test both directions. Reproduction recipe at the end of E-052
  (`freeMoves.orbitReport(..., free = freeMoves.REDUCED)`). · experimentalist ·
  open
- **T2** · Why do offsets pair by a reflection (F-053)? Each orbit is closed
  under the relation dual; find the mechanism that sends `c@o` to
  `c@(s(c) - o)`, and a formula for the centre `s(c)` in terms of the core. ·
  theorist · open
- **T3** · The exceptions to F-053 and H-020: is `3346` (five orbits at
  `n = 14`, no pairing) and `4056` (0 with 1, 2 alone) real or an artefact of
  a cap, the walk, or the gauge (F-050, F-052)? Same question for the six
  H-020 failures of E-051, all words whose slide at 13 is one or two
  offsets. · skeptic · open
- **T4** · H-020 as a theorem: the rule table acts the same on every interior
  position and only anchored and edge moves see the ends (F-051). State the
  lemma precisely and prove it, or find the rule that breaks it. · theorist ·
  open
- **T5** · H-015, the Coxeter guard: the cheap test first -- every step the
  guard admits at `n = 6` and `n = 7`, checked by a second invariant
  (see H-015 for candidates; R-008 on why Avella-Alaminos-Geiss does not
  apply directly). What does the literature offer that is computable here? ·
  scholar, skeptic · open
- **T6** · H-017, a quipu in the class carries more relations than cords:
  look for the invariant that counts relations against cords (the Euler form
  is the suggestion in H-017); test it on the recorded classes at small `n`. ·
  maverick, theorist · open
- **T7** · H-010, overlap reducible only at an end: a proof from step 7 of the
  procedure (arXiv:2112.08129, summary in `research/literature/`). The probe
  search is too long for a round; do not run `probe.py --steps 7`. · theorist ·
  open
- **T8** · Tooling for T1/T3: a command (a `batch.py` task or a small script)
  that prints E-052's report -- for a core at a length, the orbit of each
  offset under the reduced walk, its size, which offsets and which mirrors it
  holds -- with a test pinning `45` at 13. · toolsmith · open
- **T9** · H-019 and H-013 need long runs (samples at 17+, `merges.py 10
  --depths 5 6 7 8`); not for a round. A persona may write the exact overnight
  proposal, sized with `--plan`, for the human to run. · any · parked for
  rounds

Suggested first round (the chair may change it): experimentalist on T1,
skeptic on T3, scholar on T5.

## Awaiting revision

<!-- submission path · author · what the referee asked for -->
(none)

## Requests between personas

<!-- "toolsmith: experimentalist needs X" -- picked up by whoever is called next -->
(none)

## Rota

<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | - | - |
| theorist | - | - |
| skeptic | - | - |
| scholar | - | - |
| toolsmith | - | - |
| maverick | - | - |
