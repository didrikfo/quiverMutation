# Round 043 -- proceedings (ordinary)

Worked: skeptic, experimentalist, toolsmith. Referees: theorist (skeptic), skeptic (experimentalist), maverick (toolsmith). All three verdicts: minor revision. Round-042 questions were unanswered and decided by this chair (step 0.5): agenda kept, item 3 parked (no PDFs), no overnight.

## skeptic -- key-preserving J != 0 steps at n = 7 c1, c2 (E-145)
Claim: class 0 at n = 6, 7 reproduces E-143 with no c_2 != 0; at n = 7 classes 1 and 2 there are gate-admitted J != 0 steps with D = 0 (13 of 67 and 9 of 64), so E-140's "none keeps the key" fails there; the orbit relation fails on most H1/H2 steps in c1, c2. Referee (theorist): counts reproduce exactly; four class-2 steps rebuilt by hand (gate True, tiltingPlus False, Cartan incongruent, keys equal, child key computed on the child). Gaps: class-1 steps single-script; E-140 not re-run at matched budget; class numbering vs E-140 unchecked; title too broad. Decision: **accept with qualifications**. The most consequential result of the round: the "J != 0 moves the key" law (E-138, E-141) is a class-0 statement, and the key guard is no evidence for H-015 off J = 0 steps. Promoted as **E-145** with the limits written in; I did not edit E-140 or H-015 beyond this cross-reference.

## experimentalist -- n = 8 tally for E-143 (E-146)
Claim: no n = 8 J != 0 step has the E-143 shape; Q's lowest term is x^3 in some class-1 steps; E-143's "orbit relation never absent" fails at n = 7 c1. Referee (skeptic): n = 8 c1 reproduced at 150 s; counts per record, not per distinct step (25 distinct, 27 steps, 5 distinct in c1); vacuity of s, c_2 is definitional; steps not saved. Decision: **accept with corrections** (negative/null result, shallow sample). Promoted as **E-146** with distinct counts, the definitional caveat and the E-143 correction stated; "x^2 is not a law across n" softened to the explored n = 8 c1 prefix.

## toolsmith -- reverse positive control and lost edges (E-147)
Claim: the reverse search finds paths of length 2-4 (12/12 each); the 154 of 1500 lost edges are never lost to a filter but land on a different same-class algebra. Referee (maverick): diagnosis reproduced exactly; but the gate was kept in the w-search (claim "gate off" is wrong), exceptions are swallowed, and the control is shallow (depth <= 4 vs lost edges at depth 8-9; 12/12 at k = 3, 4 is improbable at a flat 10% loss), so it does not show deep reach; "closed" column wrong for one k = 4 run. Decision: **accept with corrections**; promoted as **E-147** with the gate wording weakened and the depth limitation stated. Overruled nothing.

## Questions for the steering committee
1. **Agenda:** keep the round-040 agenda, with item 1 now = why key-preserving J != 0 steps exist at n = 7 c1, c2 (which orbit data give D = 0; rebuild the class-1 steps; re-run E-140's command at 500 s) and what that means for H-015's use of the key guard. Recommend keep.
2. **Overnight:** the 3 h reverse job (`toolsmith_n6meet.py --tilting-only --reverse --hits 0,4,13`) is not added: the control is shallow only. Recommend no, until a loss-by-depth tally and a deep control exist.
3. **Literature:** still no PDFs for arXiv:1009.3370 / 2509.12983; item 3 parked. Recommend keep parked.

## Decisions taken for the steering committee
- round 042, question 1 (agenda): keep the round-040 agenda; item 1 = H1/H2 step with c_2 != 0, n = 8 tally -- decided by the chair of round 043; no answer from the human.
- round 042, question 2 (literature): no PDFs supplied; item 3 stays parked -- decided by the chair of round 043; no answer from the human.
- round 042, question 3 (overnight): none -- decided by the chair of round 043; no answer from the human.
