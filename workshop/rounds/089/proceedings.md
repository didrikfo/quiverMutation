# Round 089 -- proceedings (ordinary)

Round 088 left no questions; STEERING has no new answer. Merge of origin/main: already up to date.

## theorist -- loops of degree >= 1 at i are not reachable (T10/T5)
Claim: after one Oppermann step no loop of degree >= 1 exists, so E-235's open case never arises for the 25 failing children. Referee (skeptic): **major revision**. Classification (a)-(d) verified against the v1 source and all counts reproduced; but parents with relations carry degree-1 relation arrows (`const.usual_algebra`), and a one-step counterexample (0->1->2 zero relation, omega: 0->2 degree 1, mutate at 0) gives a degree-1 loop at vertex 2. What survives: a loop at i needs a loop in Q at i, so for parents acyclic as quivers with relations there is no loop at i. Decision: **revise** (all five required items are doable in one sitting; no response step because the verdict is major). Recorded as run E-248 with the narrowed wording; not a finding. Round 090 revision.

## maverick -- three undecided periodic groups are non-FCY on projectives (T3/T8)
Claim: R_k > n for every k <= 254, 74, 218 in 45050400, 050404400, 040404000, through the bimodule-resolution method. Referee (scholar): **accept**; reproduced the three series and a Serre check. Decision: **accept**. Promoted as E-246. Optional gaps (second prime, itemised bar agreement) go to the board.

## experimentalist -- |orbit(B_n)| = n 2^(n-5) at n = 13..16 (S-1)
Claim: fit holds, orbits closed, A_n not in the orbit (trivial). Referee (toolsmith): **minor revision** (mirror-closed column, caveat in title). Response (step 3.5) did all required items and added n = 16; I re-checked the table row n = 16 against the recurrence (32768 = 2*15360 + 2048). Decision: **accept**. Promoted as E-247.

## Consequences
- E-246 completes the rank-growth sweep of E-245: six of eight periodic indefinite groups tested (not D_10, not the 535 group) fail S^a ~ [b] on projectives, a NO beyond Cartan. No H-, F-, R- status changes; H-015 unchanged. T3/T8 next: a group-level sweep at n = 11, 12 with the cheap method, and whether growth rates split derived classes.
- E-248's counterexample means E-235's "loops of degree >= 1 at i" case is still open until the acyclic-with-relations check is done; the theorist's closure is not on the record.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- Overnight: none; the S-1 job ran in about 70 s and nothing needs more than 10 minutes per command.
- Agenda: round-088 proposed agenda unchanged; item 1 gains the theorist's revision, item 2 gains the sweep of all periodic groups with `maverick_bimod.py`.
