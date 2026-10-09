# Round 038 -- proceedings (ordinary)

Worked: experimentalist, toolsmith, scholar. Referees: theorist (experimentalist), skeptic (toolsmith), experimentalist (scholar). Step 0.5: round 037's questions were unanswered; decided and recorded in `STEERING.md` (keep the round-036 agenda; no overnight).

## experimentalist -- d >= 3 with J != 0 on the n = 8 c0 walk
Claim: 9 rows / 5 algebras, all out(i) = 3, two new (3,1) rows without parallel arrows; out(i) = 3 necessary in the sample, not sufficient. Referee (theorist): minor revision; every d >= 3 cell reproduces, d = 2 counts move with the time cap; asks for per-row listing, conditional ratio, cwd note, isomorphism caveat. **Decision: accept with qualifications** (all five asks adopted as Limits in the entry). Promoted as **E-137**.

## toolsmith -- the 16 fans meet an LNA; maverick_single fix
Claim: all 16 E-134 fans share 25 algebras with the LNA-side forward graph; n = 6 classes do not close; `maverick_single` label bug fixed. Referee (skeptic): minor revision; reproduces (27 shared, cap drift); None-artefact ruled out; equivalence rests on the guard (H-015), not the gate; no `tiltingPlus` check of any step; E-087 caveat untested. **Decision: accept with qualifications**, wording corrected to the guard. Promoted as **E-136**, with a note on E-134 that its "no LNA reached" is superseded. Test `tests/test_toolsmith_single.py` added. Consequence: E-131's d = 2 is a forward-walk statement, not a derived-class one.

## scholar -- Cartan defect formula, no known invariant bounds d_i
Claim: `C_B = r C_A r^T + H`, `H_{vi} = dim J_i`; d_i is a Cartan entry, not a class invariant. Referee (experimentalist): minor revision; reproduces, 0 failures in 28 000 walk steps (n = 6, 7), but largely recorded as E-095 / E-129, H from the repo's own `perI`. **Decision: accept with qualifications; "new" softened to "new as a closed formula with derivation"**. Promoted as **E-138**. No hypothesis status changed.

## Questions for the steering committee
1. **Agenda:** keep the round-036 agenda. Item 1 is reshaped again: since the 16 fans (d >= 3, J != 0) lie in the LNA class at n = 6 (E-136), ask for an independent check of one meeting path with `tiltingPlus`/Ladkani 2.3(c) (skeptic), and why out(i) = 3 (theorist). Recommend keep.
2. **Overnight:** `toolsmith_n6close.py --budget-hours 3 --only 0` was proposed by the toolsmith (does class 0 stop growing?). Recommend not yet; an `isTilting`/`tiltingPlus` check of the meeting path comes first.

## Decisions taken for the steering committee
- round 037, question 1 (agenda): keep the round-036 agenda -- decided by the chair of round 038; no answer from the human.
- round 037, question 2 (overnight): none -- decided by the chair of round 038; no answer from the human.
