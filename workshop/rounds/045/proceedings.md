# Round 045 -- proceedings (ordinary)

Worked: experimentalist, toolsmith (T10 guard audit, a and b), theorist (T7, breadth slot). Referees: skeptic (experimentalist, toolsmith), maverick (theorist). Step 0.5, STATE rewrite and the H-015 ledger pass followed the special request in STEERING (retrospective 1).

## Submissions
- **experimentalist** (`experimentalist.md`): the 10 F-041 n = 8 merges also meet in a tilting-only walk; all 2396 gate-admitted edges are J = 0 with the key kept; guarded depth-4 LNA-started edges at n = 6..8 are all J = 0 (n = 8 partial). Referee (skeptic): minor revision; reproduced n = 8 depth 3 and n = 6, library and BFS agreed on 2 pairs, "not new information, follows from E-084". Response (step 3.5): all four items done, library `meetingPoints` agrees on 10 of 10. **Decision: accept, narrowed** (scope: depth <= 4, forward from LNAs, F-041 pairs only; tilting-only meeting is implied by the edge identity). Promoted as **E-148**.
- **toolsmith** (`toolsmith.md`): an added J = 0 / `tiltingPlus` check costs 8-11% per step; key-keeping steps failing `tiltingPlus`: 16 of 80 978 (n = 7 c1), 9 of 79 143 (c2), at parent depth 7-8, none at n = 8 samples. Referee (skeptic): minor revision; reproduced the n = 7 c2 counts and the c1 cost row; the E-084 conflict is not a contradiction (E-084 stopped at its first rejecting level); "docstring wrong" rested on the untested premise that the children leave the class. Response: five of eight items done (E-084 reconciled, depths, taint tally, claim restated, ratios only); not done: n = 8 reruns, E-084 class-index mapping, out-of-class test. **Decision: accept, narrowed** to what the done items support. Promoted as **E-149**.
- **theorist** (`theorist.md`): no proof of H-010; a "run-of-three" lemma (a bystander lowers an interior heavy pair only if it shares >= 2 arrows); recommends closing T7. Referee (maverick): **major revision**; k = 2 table reproduces, no counterexample, but the lemma largely restates F-022's run-of-three table, "22 of 22" is 21, "holds at 3 mutations" is unsupported (k = 3 hit the cap), the closure is unsupported. **Decision: revise** (the referee is right; the closure recommendation does not follow). Nothing promoted. Onto *Awaiting revision*: items 1-7 of `theorist.review.md`, in particular `(8:3)(9:3)(10:4)` at k = 3 and the left bystander.

## Consequences
- **H-015** moves from SUPPORTED to **OPEN** (ledger pass, `research/HYPOTHESES.md`): its sufficiency reading is refuted for *tilting* at n = 7 c1, c2 by E-145 and E-149; derived-equivalence of the failing children is untested, so no R-entry. The status line is one sentence; history moved into the body.
- **`search.mutationSearchDepthFirst` docstring** (`quivermutation/search.py`, line 322: "`coxeterGuard` is what makes the walk a walk in one derived class") is not supported by E-149 as a statement about tilting, and unproven as one about the class. Not edited this round: the rewording waits on the skeptic's out-of-class test (thread T10 item). Thread T10 stays open for it.
- E-084's "0 of ~1.3e6" is bounded to the walk depth it covered (stopping level 5-7); it does not extend to depth 7-8 (E-149).
- E-140's "none keeps the key" remains superseded for n = 7 c1, c2 (E-145).

## Decisions taken for the steering committee
- round 044, question 1 (agenda): proposed round-044 agenda approved, T10 ahead of it -- decided by the chair of round 045; no answer from the human.
- round 044, question 2 (overnight): none -- decided by the chair of round 045; no answer from the human.
- round 044, question 3 (literature): parked; arxiv.org is denied by the network policy -- decided by the chair of round 045; no answer from the human.

## Questions for the steering committee
1. **Promote `tiltingPlus` into `quivermutation`** (and add `tiltingGuard = False` to `mutationSearchDepthFirst`, cost 8-11% per step)? Only the human can overturn STEERING q3 of round 001. Recommend: yes, as an opt-in keyword, default off, once the skeptic's out-of-class test of the E-145 children is in; until then no.
2. **Allow `arxiv.org` in the environment's network policy or supply PDFs of 1009.3370 and 2509.12983** (only the human can). Recommend: yes.
