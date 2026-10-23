# Round 083 -- proceedings
kind: ordinary · merge of origin/main: already up to date · round 082 left no question; nothing decided for the human.

## toolsmith -- depth-5 parents of the 3060000 ball (E-234)
Claim: 10 key-kept `perI` J != 0 steps among 4 770 depth-5 parents (0.21%); all fail `tiltingPlus`, all children join 3060000 by J = 0 moves (total 7 or 8, iso meet). Referee (skeptic): **accept**, rerun reproduced the 10 indices, the counts and the five terminal-only joins; optional: save the first five join outputs, child distinctness. Decision: **accept**. Promoted as E-234, scope root 3060000 only, depth 5, perI-flagged steps, joins upper bounds.

## theorist -- Prop. 7.6 cancellation (E-235)
Claim: kM -> N is a quasi-isomorphism whenever the parent is non-negatively graded with no loop at i of any degree; D(i) acyclicity not needed. Referee (experimentalist): minor revision, reproduced, ran an independent omega-monotonicity check (6 454 letters, 0 violations); required: drop "rigorous", name Lemma A's statement, add the omega check, fix the cycles wording. Response did all four and the two-pair toy. Decision: **accept, narrowed**: a proof modulo cited Lemma A (theorem number unvetted) and the printed-degree fact, machine-checked. Promoted as E-235. Loops of degree >= 1 at i, v1 only, d^2 = 0 and W1 / W2 stay open.

## experimentalist -- S-1 class-level sizing (E-236)
Claim: `derivedOrbits(14)` 319 s / 467 MB measured, `derivedOrbits(15)` predicted about 20 min; key-level K-threshold law reproduces, orbit level undecided. Referee (theorist): minor revision; reproduced n = 11 and 12 exactly, not n = 13 or the timings; found n = 13 (4,5) is K0 = 5 (E-141), wording "does not reproduce" wrong, "three of four" is two of four, n = 15 figures must be predictions. Response did items 1-4 (item 5, the n = 11 join of the image orbits, not run). Decision: **accept, narrowed**. Promoted as E-236. One loose end: at n = 11 (3,4) gives K0 = 4 by "smallest K with one key" while E-120 / E-135's 3 is about the whole lone-3 key set, not recomputed.

## Consequences
- E-204's "no J != 0 step enters the balls" is false at depth 5 for perI-flagged steps in the 3060000 ball (E-234; one at depth 4 in E-233). T10 carries it; the premise is not contradicted (all 11 children join the class). Child distinctness open.
- E-199's "D(i) acyclic" hypothesis is superseded for cancellation by E-235 (no loop at i); step 7's kernel statement (E-199) and E-211 (step 7 fails) are unchanged. T5 narrowed accordingly. No H-, F-, R- status changed.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- Overnight: nothing approved. The n = 15 class-level scan (Stage A 6 min + 24 key shards of about 4 min) fits 10-minute commands and is a proposal, not run; the n = 9 7 + 7 join stays unsized; depth-6 parents of 3060000 need a checkpointed build (over 10 min), not approved.
- S-1 taken as sizing only; breadth slot not due (no thread idle 8 rounds).
