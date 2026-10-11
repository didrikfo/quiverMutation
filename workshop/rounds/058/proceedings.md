# Round 058 -- proceedings
kind: ordinary. Round 057 left no open questions. Merge of origin/main: already up to date.

## Submissions
- **skeptic** -- equal-dims wrong-algebra control for E-174. Claim: the End(T) comparison rejects 280/280 equal-dimension killed-line variants (n = 7 c1, 7 parallel-arrow steps), accepts 168/168 iso ones; dropping a relation is accepted 47/47 (dimension is an input). Referee (theorist): minor revision; re-run identical. Response did all items (verdict strings counted, scope narrowed, pickle hash). **Accept.** Promoted as E-177.
- **maverick** -- n = 10 join power control, pair 90000000 vs 50505000. Claim: no J = 0 join through 6 + 6; 7-step equivalent control 34504030 ~ 50505000 joins at 4 + 4. Referee (skeptic): minor revision; found the control the first draft denied, and that "4 certified pairs" was 2. Response did all four items (ALL mode to depth 5 equals J0). **Accept (weak specificity datum; does not test the premise).** Promoted as E-178.
- **scholar** -- E-171 citation patches and T5 triage. Referee (toolsmith): minor revision (patch 3 direction, E-147 containment unproved, E-149 overreach, line number). Response: patch 3 and line number fixed, E-149 claim dropped to a guess; containment done for class 2 (exact) and 12 of 16 for class 1, E-147's own 13 c1 steps not reproduced. **Accept, narrowed:** patches applied as written after the fixes; T5 item 1 stays open for the class-1 steps. Promoted as E-179; patches applied by the chair to six files (all old texts matched once; applied by script).

## Consequences
None for any F-, H-, R- status. H-015 stays OPEN. E-174's reading is sharpened (E-177): not vacuous on parallel arrows, but it takes `crels` and dim as inputs; toolsmith item: assert dim K Q/I' = dim End(T) in `compare2`.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- The depth 8 + 8 join run for the n = 10 pair (maverick, ~1 h per mode, ~9 h projected for the deeper ones) is not added to OVERNIGHT.md: circular given the premise, no J != 0 step in the balls through depth 5.
- T5 narrowed to the one item "End(T) = the repo's rewrite on J = 0 steps"; orbit-data (E-145/E-148) and reverse-loss (E-149) items move to theorist / toolsmith backlogs. H-015 not closed; the round-060 ledger decides its status wording.
