# Round 067 -- proceedings
kind: ordinary. Called toolsmith (T10/T5), theorist (T5), experimentalist (S-1 item 2 = agenda item 3). Referees: skeptic -> theorist, theorist -> experimentalist, maverick -> toolsmith. All three: minor revision, all answered in step 3.5 (checked against the reviews; the toolsmith's response ran the full suite and the end check).

## toolsmith -- every printed-path parent acyclic; cyclic-input warning
Claim: no loop or oriented cycle among 405 mutated algebras (n = 7: 370, from E-168, not rerun; n = 10: 35) and 49 ends; `CyclicQuiverWarning` added. Referee: minor revision (end check, non-vacuity control, dual-invariance sentence, full suite). Decision: **accept**; response did all four. n = 7 half is a rediscovery of E-168, said so. Promoted: **E-198**; GLOSSARY `CyclicQuiverWarning`. Library change kept (warns only, values unchanged; full suite 1426 passed).

## theorist -- g = gd and Lemma L without acyclicity; D(i) criterion
Claim: proofs (A), (B), (C); repo step 7 equals Oppermann's kernel whenever D(i) is acyclic; 7 389 pairs with a cycle upstream, 0 differences. Referee (skeptic): proofs check line by line, re-run identical; fixes: 'on' count 5 840, dead `thruI` counter, headline should say the upstream test checks code only. Decision: **accept, narrowed** as the response did. Promoted: **E-199**; GLOSSARY D(i); H-015 status line pointer. Does not prove J = 0.

## experimentalist -- n = 9 -> 8 deletion table
Claim: none of class / i / left-right coverage fixes the image class beyond E-114's free-stretch rule; key = class at n = 8. Referee (theorist): all counts reproduce; mostly E-114 / E-146 / E-194; wording past the four predictors; class letters unstable. Decision: **note** (rediscovery plus a table); the response narrowed the title, added cell counts (two classes rest on 4 and 2 cells), made letters deterministic, leave-one-out 62.6%. Promoted: **E-200** (run record).

## Consequences
None for F-/R- statuses. H-015 stays OPEN. The agenda item "every parent on a printed path is acyclic" is closed for these paths; the acyclic-parent restriction of G2 (E-193) is replaced by "no oriented cycle meets D(i)" (thread T5 updated). Open: cyclic children as later parents off the printed paths; a cycle meeting D(i).

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- No unanswered questions from round 066 (it left none). No overnight run proposed.
- Round 068 is a conference (multiple of conference_every = 4).
- S-1 item 2 taken at n = 9 -> 8 (not taken in 066); no breadth slot (no thread idle 8+ rounds; T6 oldest at 6).
