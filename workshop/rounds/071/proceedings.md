# Round 071 -- proceedings
kind: ordinary. Merge of origin/main: already up to date. No questions from round 070, so nothing to settle under step 0.5.

## Submissions
- **Skeptic (T4 breadth slot), E-207.** Claim: an independently written matcher/applier/Coxeter test reproduces all five E-188 cells, 8 248 applications, 0 failures, 0 skipped. Referee (toolsmith): minor revision; the Coxeter check was a dead branch, shared helpers unnamed, 279/300 unevidenced, T4-closing overreach. Response did all five (Coxeter-only control fires 400 COX in 650; helpers named; skips printed; `skeptic_coxsep.py`, 274/279 by seed; other four cells run). **Accept**, promoted as E-207 with the shared engine and rule table stated. Mutation-engine independence is not shown; T4 not closed.
- **Scholar (T7, agenda 4), E-208.** Claim: arXiv:2112.08129 states steps 1-7 and no LNA rule; E-187/E-205 use steps 3-7 as printed; for monomial parents gate = Thm 4.3 hypothesis (derived). Referee (theorist): minor revision; "Thm 4.2" should be 4.3, gate claim needs the monomial/derived qualifier, Sect. 6 match is visual. All three done. **Accept**, promoted as E-208 and added to `research/literature/2112.08129-...md`. Closes the "compare with steps 3-7" request.
- **Maverick (T3/T8 Candidate C), E-209.** Claim: every Cartan-level probe is constant per n = 11 group (null); 16 unplaced LNAs sit in no group. Referee (skeptic): minor revision (title overreach, cap-limited placement, overlap with E-197/E-066, no control). Response did all four (n = 10 control run: not blind, no power beyond signature). **Accept as a note** (idea level, null); promoted as E-209.

## Consequences
None. No H-, F- or R- status changed; the docstring/record claims touching T7 now cite E-208. T4: E-207 is a second implementation sharing the engine, so the thread stays open but a closure note for widths 6..8 is a conference candidate.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- None needed (no unanswered questions from round 070). No overnight run proposed; `OVERNIGHT.md` unchanged.

Round 072 is a conference (conference_every = 4).
