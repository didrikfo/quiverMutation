# Round 095 -- proceedings
Kind: ordinary. Merge of origin/main: already up to date. Round 094 left no questions; STEERING has no new answer (step 0.5: nothing to decide).

## Submissions
**toolsmith** (rediscovery; T10/T5, agenda item 1). Claim: the 2400230(8) vs 3304000(6) union join was already run in round 093 as E-255 (MISS, move-set statement, no bearing on J = 0, since E-238 certifies the pair equivalent); only `--plan` was run, after stating what a MISS means. Referee skeptic: minor revision (true, not new; 3 small items). Response: all done (8427/8063 isomorphism output located at `workshop/rounds/093/toolsmith_extra_out.txt` line 5; claim now opens with "no new computation"; E-244 annotation requested). Decision: **note**. Nothing new to promote; I annotated the "8427 not run" line in E-254's scope with a pointer to E-255. Agenda item 1 is closed as done by E-255; no larger radius (a MISS cannot become a NO by repetition).
**theorist** (proof, partial + negative; T4 breadth slot). Claim: H1 of H-020 split by proof coverage: 110 of 414 floating rules (class A) replay exactly under E-187/E-205 one-step rules; the other 304 (class B) leave the LNA class at the first mutation, where E-187 is silent. Referee skeptic: minor revision (reproduced 9 396 cases / 0 violations, A = 110, B = 304; "follows by locality" overstated; label-swap empirical). Response: retitled, wording narrowed, widths 9 and 10 at w+4 pass, A/B split unchanged at w+4; width 11 at w+4 not done (`[next round]`, 2 rules). Decision: **accept, narrowed** -> E-262. H1's status unchanged.
**scholar** (negative, search-bounded + citation audit; T3/T8 and T6 literature). Claim: no source found defining H-017's "relations minus cords"; Tits form = Euler form as matrices iff gl.dim <= 2; citations of E-251/E-253/E-257 match their sources except the rank-1 step of twisted FCY; Serre-dimension slope re-reads E-253. Referee maverick: minor revision (parity misstated as (-1)^b, should be (-1)^(a+b); b/a comparison already closed by E-253). Response: all five items done. Decision: **accept, narrowed** -> E-261 (the Tits-form sentence is UNVERIFIED against Bongartz; the Serre-dimension row is not a new test).

## Consequences
None. No H-, F-, R- status changes; no library code or docstring contradicted. Open points carried: width 11 at w+4 (T4); twisted-FCY rank-1 step at Herschend-Iyama's definition (skeptic, T3/T8).

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- No overnight run approved (none proposed that is needed; toolsmith sized the union join at about 20 min in 10-minute commands).
- Agenda kept. Item 1 closed (done by E-255); T10 now waits for a non-move-set certificate. S-1 not taken (last 089; no slot spare). Round 096 is a conference.
