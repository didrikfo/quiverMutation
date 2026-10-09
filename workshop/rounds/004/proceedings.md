# Round 004 -- proceedings

Step 0.5: round 003's one question (add the n = 12 and n = 14 censuses to `OVERNIGHT.md`) had no answer; decided by the chair as recommended (yes), recorded in `STEERING.md`, entry added to Menu 4 (test_overnight_doc: 110 pass). `main` merge: no-op. Note: round 004 is a multiple of `conference_every`, but round 003's STATE said `ordinary` and I followed it; round 005 is a conference.

## theorist -- T2
* **Claim.** The interior/end-touch split is a committed column (13/13, 17/21 as E-062) and explains none of the 7 failures; for `33x`, `k = 2x`, `d = x - 3`, n = 13..17, x = 3..6.
* **Referee (skeptic): minor revision**; all three scripts byte-identical. Wanted: `337` run, number of cores with a unique centre, the parity reading of `4046 5046 5056` marked unverified.
* **Decision: accept.** I ran `337` at n = 15, 16, 17 (`k = 14`, `d = 4` each; low power at 15), counted unique centres (59 of 109) and marked the parity reading unverified in E-063.
* **Promoted.** E-063; H-021 status line.

## experimentalist -- T1
* **Claim.** The 12 cores of E-062 keep `k`, `d` at n = 15 (12/12 fit) and 16 (9/12); the three others have an unmerged equal-size middle pair; parity does not carry over.
* **Referee (scholar): minor revision**; re-ran `46` at 16 and the table, identical. Wanted: cite E-060 by identifier, say "no fit" means unmerged under the round-002 rule, note the size 20300 recurs from E-060.
* **Decision: accept**, all three applied in E-064 (the 20300 question is stated as not checked).
* **Promoted.** E-064; H-021 status line.

## maverick -- T6 (first submission)
* **Claim.** H-017 not refuted at n = 9 to depth 6, but the Coxeter polynomial and Euler form do not predict (cords, relations); the Euler-form signature separates "outside every quipu class" for n = 8..11.
* **Referee (skeptic): minor revision**; reproduced n <= 10 signatures, the n = 13 bad quipu, `c_{n-1} = 1` (plus 3000 random monomial trees). Not reproduced: n = 12 (over 10 min). Trace fact for hereditary trees is Happel's (`literature/2509.02375`); the criterion is verified for n = 8..11 only; the n = 13 "failure" is a proof-method break.
* **Decision: accept with those wordings** (all four required items done in E-065 as wording or marked unrun; I did not run n = 12). Direct search from 16 candidates recorded as not evidence.
* **Promoted.** E-065; H-017 status line.

Glossary: nothing new.

## Decisions taken for the steering committee
- Round 003 question 1: n = 12 and n = 14 censuses added to `OVERNIGHT.md` Menu 4 (recommended option; no answer from the human).

## Questions for the steering committee
1. Maverick proposes two H-017 runs for `OVERNIGHT.md`: `maverick_reached.py 9 7` on `3033030` and `4444400`, and depth 4 over the 262 outside LNAs at n = 10. The n = 10 run cannot decide anything until the search has a positive control (a certified member found from its LNA at its path depth). **Recommend: approve the n = 9 depth-7 run only; toolsmith builds the control first.** If unanswered the next chair takes this.
2. Round 005 is a conference (STATE had wrongly said `ordinary` for 004). Recommend: keep it a conference.
