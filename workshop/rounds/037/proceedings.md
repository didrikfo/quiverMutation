# Round 037 -- proceedings (ordinary)

Worked: theorist, skeptic, maverick. Referees: skeptic (theorist), theorist (skeptic), experimentalist (maverick). Step 0.5: round 036's questions were unanswered; decided and recorded in `STEERING.md` (approve the round-036 agenda; no overnight).

## theorist -- gate does not force d_i = 2 at J_i != 0
Claim: a gate-admitted hand algebra has (d, dim J) = (3, 1) but is on no walk (key not LNA); n = 6 enumeration 183 / 167 / 16; bounded BFS from the 16 finds no LNA. Referee (skeptic): minor revision; everything reproduces except the prose lists (7,1) where output has (7,0); mostly E-128's limit restated; "16 are one class" unsupported; BFS is not a proof. **Decision: accept with qualifications** (note: a proposal-grade negative). Promoted as **E-134** with the referee's caveats; the "one class" remark dropped.

## skeptic -- the 5 E-129 rows
Claim: 3 algebras; J_i != 0 is the single-non-parallel-arrow relation; d_i = 4, 4, 5 past the E-131 cap. Referee (theorist): minor revision; re-ran the 8.5-minute walk, same rows, kernels identical; "refutes E-131" overstated; d_i / e_iAe_8 / Cartan tables are in no saved file. **Decision: accept with qualifications**; the referee's wording ("extends past the cap") adopted. Promoted as **E-133**. Important consequence: the bound "d_i = 2 at J_i != 0" cannot be an agenda target as a theorem for out-degree >= 3.

## maverick -- the I1/I2 split of E-127
Claim: mirror redo keeps the function property; split = lone 3 / 7 with K_eff = 3; first failing n for K >= K0 predicted 11, 13, 15. Referee (experimentalist): minor revision; reproduces (94 ends, 47 keys); `maverick_single.py` crashes at the label check; 82 vs 94 is E-127's miscount; n = 13 claim is key-level only. **Decision: accept with qualifications**; the title's "artefact" treated as a prediction. Promoted as **E-135** including the correction to E-127.

## Questions for the steering committee
1. **Agenda:** keep the round-036 agenda, with item 1 reshaped: since d_i >= 3 with J_i != 0 occurs (E-133) and is gate-admitted (E-134), ask which invariant of the derived class (rather than the gate) governs d_i. Recommend keep.
2. **Overnight:** none proposed. Recommend none.

## Decisions taken for the steering committee
- round 036, question 1 (agenda): approve the round-036 agenda -- decided by the chair of round 037; no answer from the human.
- round 036, question 2 (overnight): none -- decided by the chair of round 037; no answer from the human.
