# Round 086 -- proceedings
kind: ordinary. Merge of origin/main: already up to date. Round 085 left no questions; nothing to settle.

## Theorist (T5) -- E-240
Claim: Anick's Lemma A in the monomial case (X1, Y disjoint, D = 0 on arrows) has a direct skeleton-decomposition proof; 13 748 cells, 0 mismatches; shared-arrow control fails 164 of 215. Referee skeptic: **accept**; re-ran the final script, byte-identical output; proof checked, no gap. Decision: **accept**. Promoted as E-240. Caveat kept: that E-235 step 4 lands in these hypotheses is not shown, so the E-235 annotation is conditional.

## Experimentalist (T7) -- E-241
Claim: pair (8:3)(9:3), k = 4, 5 of 5 m = 2 two-bystander cells stay at overlap 2; 3 of 3 completed controls drop to 0. Referee theorist: **minor revision** (joint run of three per cell; "survives" -> "not contradicted"). Response done: every cell has a joint run of three because the pair is a run of four; in 2 of 5 the bystanders alone also form one and still stay at 2; claim narrowed. Decision: **accept, narrowed**. Promoted as E-241. Open: k >= 5, m >= 3, other shapes, (10,4)(5,2) control (timeout).

## Maverick (T3/T8) -- E-242
Claim: n = 12, 8 key groups with indefinite periodic LNAs (8 745), all 8 737 non-representatives Cartan-congruent to the representative; no real NO control exists at n = 12. Referee skeptic: **minor revision** (count 8745, "representative" among indefinite periodic members, P checked on (G, Phi)). Reproduced exactly (445 s, 8 737 YES). Response done; pickle path made reproducible (`$MAVERICK_PKL`). Decision: **accept** as a null result. Promoted as E-242.

## Consequences
None for any H-, F-, R- status. E-235's "modulo the cited Anick lemma" is discharged for the monomial step only conditionally on the step-4 hypotheses check (thread T5). The Cartan-level NO still lacks a key-equal incongruent control through n = 12 (T3/T8).

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- Overnight: none approved. Proposed and left unapproved: n = 13 Candidate C run (~40 min, maverick); 12 remaining k = 4 m = 2 cells plus (10,4)(5,2) (~4 h core time, experimentalist).
- Agenda unchanged.
