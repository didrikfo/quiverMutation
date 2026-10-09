# Round 018 -- proceedings

Ordinary. Worked: skeptic, scholar, maverick. Referees: theorist (skeptic), experimentalist (scholar), skeptic (maverick). All three came back **minor revision**; the points were scope and wording that I applied while promoting (no new runs needed), so all three are accepted with the referees' limits written into the entries. Step 0.5: STEERING answers for round 017 recorded (below).

## skeptic -- row sets of the `444` orbit at n = 12..17
- **Claim:** among 4-letter words (letters <= 9, a 4, >= 4 placements) 6, 9, 12, 15, 18, 18 are split across the `444` orbit at n = 12..17, each with exactly one placement in it; `3334`, `2455` have none; E-086's "0 partial" holds only for merged words, where it is true by definition.
- **Referee (theorist):** reproduced byte-identically (n = 14, 16) and recounted all six n. Errors: 55 walked words not 56, "bounds the sum" explanation of the 120 plateau, quantifier too wide, n = 16/17 lists are one family.
- **Decision: accept**, with the quantifier narrowed and the errors corrected. Promoted as **E-091**; H-021 status line updated; glossary: *split word*.

## scholar -- Cartan congruence vs `tiltingPlus`
- **Claim:** on all gate-admitted steps tested (n = 5..7, 807 non-tilting steps) Cartan congruence fails exactly where `tiltingPlus` fails, the discrepancy being row k, off-diagonal, minus the kernel dimension of `p -> (p beta)_beta`; so congruence is `tiltingPlus` read through the rewrite, not a second criterion.
- **Referee (experimentalist):** E-078 and n = 5 reproduce exactly; the n = 7 wall-clock-capped run gives different counts (24 187 / 44 241 / 120) with the same pattern. No histogram of dim ker or count of distinct parents; derivation covers one direction at (k, i) only; mechanism read off the data.
- **Decision: accept as an observation, not a theorem** (title softened in the entry; the limits list the missing histogram and the cap dependence). Promoted as **E-093**; H-015 status line updated. I did not run the histogram; it is a request.

## maverick -- positive control for `MONO=1`; producing relations
- **Claim:** all 2376 n = 8 cord members (six LNAs, length <= 5) carry a sum relation, none is monomial; `MONO=1` finds no monomial cord at n = 4 (depth 8) and n = 5 (depth 7); the filter is positive on a hand-built seed; no positive control at n = 8.
- **Referee (skeptic):** numbers reproduce; the n <= 5 negatives are informative for only 6 of 14 (n = 5) and 1 of 5 (n = 4) LNAs; n = 8 is depth 5, shallower than E-087's 6; "covered" is a union of supports.
- **Decision: accept as a scoped negative**; the entry states all three limits. Promoted as **E-092**; H-017 status line updated.

## Questions for the steering committee
1. **Overnight:** none proposed. Recommend: no. (The depth-8 n = 8 class 2 walk still wants a checkpoint in `scholar_walk.py`, toolsmith.)
2. **Agenda:** approve the round-016 agenda in `STEERING.md`? Recommend: yes, unchanged. Round 019 is ordinary; 020 is the next conference.

## Decisions taken for the steering committee
- Round 017, question 1 (overnight): none; the depth-8 n = 8 class 2 walk needs a checkpoint first -- decided by the chair of round 018; no answer from the human.
- Round 017, question 2 (agenda): approve the round-016 agenda unchanged -- decided by the chair of round 018; no answer from the human.
