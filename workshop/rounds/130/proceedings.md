# Round 130 -- proceedings
Cycle 2, phase 3 (round 2 of 4). Ordinary. Step 0.5: round-129 question 1 decided by the chair (see below).

## Submissions
- **toolsmith** (T10, F-041 with `refuseNonTilting`): at n = 10, 11, depth 3 + 3, the option reaches the same quivers as the default in all 3 634 balls; no join or class merge lost (206/304, 924/1 513 joined). Referee maverick: minor revision (calls not counted; 12 direct rows not committed; the control fires on a feature the balls lack). Response did all three; thinness stated (0 of 172 047 counted steps has dim e_i A e_v >= 2). **Accept**, promoted as **E-337**, F-041 date line updated. The null means "the balls are thin", not that the option was tested where it could fire; n = 11 calls counted on 12 pairs only. Library change: `refuseNonTilting` passed through `quiversReachedFrom`/`meetingPoints` (default off), test added.
- **skeptic** (T4, H-029 widths 8-9): one delta per rule, 106 of 106 rules, 2 332 matches; shuffled-label control has power; the registered cross-width falsifier passes without power. Referee scholar: minor revision (falsifier wording, scope, near-s0 audit, planted control). Response did all four. **Accept**, **E-335**, **H-029 -> SUPPORTED** at the narrowed scope. T4 closes as a rule-table thread (proposal accepted); relation counts are the engine's only.
- **experimentalist** (T7, interior gap-2 control): (8,5)(10,4) with bystander (13,2) or (12,3) lowers overlap 3 -> 2 under E-330's harness, so E-330's null for that shape has power; (14,2) does not; no control for (8,4)(10,3) (0 of 91 + 1242 non-vacuous screened cells; screen covers double mutations only). Referee theorist: minor revision (claim too broad, vacuous counts). Response did all five except the capped k = 5 cell (overnight-sized). **Accept, narrowed**, **E-336**; E-330 amended: one of its three non-vacuous cells now has power.

## Consequences
- E-330's "no interior control" is partly false: a control exists for (8,5)(10,4). H-010 stays SUPPORTED (two of three cells still have no power). Thread T7 updated.
- H-029 status OPEN -> SUPPORTED; thread T4 closed.
- Round 129's question 1 resolved at depth 3: flipping `refuseNonTilting` costs no F-041 join at n = 8-11; the only known cost remains E-316's mixed F/R joins. Docstring of `mutationSearchDepthFirst` still false for the shipped walk (backlog).

## Questions for the steering committee
None. (Default stays off; revisit at the conference ledger.)

## Decisions taken for the steering committee
- Round 129, question 1 (`refuseNonTilting` default): keep it off, document the crossing, toolsmith run F-041 at n = 10, 11 first -- decided by the chair of round 130 (recommended option); done, E-337. Recorded in STEERING.md.
- No overnight run proposed. The capped k = 5 cell of E-336 and a depth-4 refusal count on the unjoined F-041 pairs are sized as overnight candidates; not requested.
