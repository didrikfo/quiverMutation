# Round 087 -- proceedings
kind: ordinary

## Settled from round 086
Round 086 left no questions (*Questions for the steering committee*: none); STEERING.md has no new answers. Nothing to apply. Merge of origin/main: already up to date.

## Submissions

### skeptic -- T4 breadth slot (rule: idle since 079)
Claim: for the stored `VERIFIED_MOVES` width-5 sequences at length 10, an independent silting-mutation engine gives End(T) = the predicted LNA in 1 008 of 1 008 applications; three controls (wrong result, wrong sequence, same-Coxeter-polynomial LNA) reject all 1 008; extra cells w5 L11, w4 L10, w4 L11 clean. Referee (theorist): minor revision (k-range, generation, saved outputs, wording of what is audited). Response done, all five items. Decision: **accept**. Promoted as **E-243**, scope "sequence-to-End audit, not the rule table". No H- status change.

### scholar -- T3/T8 Candidate C
Claim: rank of nu^k P_l: returns to n at k = 15 with nu^15 P_l = P_l[16] in the 535 group (rad^4, rad^5, 033450000); grows for every k tested in 34504030, 50505000, 334504030, so periodic Phi does not give FCY; controls pass. Referee (skeptic): minor revision; reproduced rad^4, 34504030 growth, 30302, plus an independent Serre-duality cohomology test (0 mismatches). Response: independent check adopted, undecided groups (45050400, 050404400, 040404000) removed from the claim, "linear growth" dropped, misprint dropped (already recorded in literature/1310.1557 lines 35-42), certificates saved (rad^5 and 033450000 outputs are from the original runs, not re-run). Decision: **accept, narrowed**. Promoted as **E-245**.

### toolsmith -- T10 3304000 depth 6
Claim: depth-6 scan fits in 10-minute commands; 6 key-kept perI J != 0 steps among 11 355 parents, all real, 3 iso classes, five of six children join by J = 0. Referee (experimentalist): minor revision; Hom tests, distinctness, sub-range scans and three joins reproduce; build and two big shards not independently rebuilt. Response: RANGE lines added, join lists replaced, title narrowed, 8196 and 11123 joined. 8427 not joined, build not independently reproduced. Decision: **accept**. Promoted as **E-244**; no overnight run needed for this root.

## Consequences
- E-245 contradicts the reading of "periodic Phi" as support for Candidate C's S^a ~ [b] in the three (2,0,*) groups (E-231 listed FCY status for rad^4/rad^5 as UNVERIFIED; it is now FCY on projectives, 16/15). T3/T8 updated on the board: the rank-growth invariant is a NO with power beyond Cartan data; next is to use it on the n = 11 groups.
- E-244 repeats E-234/E-239: "no J != 0 step enters the balls" (E-204) is false for both roots at depth 6 too. Already on T10.
- No H-, F-, R- status changed; no library or docstring claim touched.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- Overnight: none approved; nothing needs more than 10 minutes per command. 3060000 depth 6 sized by the toolsmith (about 12 min build in two commands, 7 min scan): not an overnight job.
- Agenda unchanged; round 088 is a conference.
