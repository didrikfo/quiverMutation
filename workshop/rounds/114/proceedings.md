# Round 114 -- proceedings (ordinary, cycle 1 phase 4: research)
Step 0.5: round 113's only question (PDFs of Nikulin 1979, Happel-Seidel, BGS, Anick 1982) unanswered; decided: citations stay UNVERIFIED, no action (STEERING). Merge of origin/main: up to date. Two pushes were refused once and went through on retry.

## Submissions
- **Toolsmith (T10, step 8)** -- Hom-verified 12-edge chain `R4 F1 R7 F1 F4 F3` / `F1 F3 R7 R4 R3 F5` joins parent and child of E-201's step 8; with E-297 all 13 steps joined (YES side, meets iso up to coefficients, generation assumed). Referee skeptic: minor revision (13/13 rests on E-297; controls test the meet iso more than Hom; replay time; A-side control), answered (A-side sample 34 of 84 rejected, 50 left). I re-ran nothing; the referee's replay reproduced YES. Decision: **accept** (narrowed wording already in the title). Promoted as **E-300**.
- **Experimentalist (S-1)** -- direct tail deletion equals the head tables at n = 17, 19, 21; off-centre lone-3 classes have K0 = min(h, n-4-h) (7 classes, a fit); E-263's "messy" (4,7) class is a digit-join writer bug. Referee theorist: minor revision, answered. I checked the writer line (`experimentalist_scan15.py` l. 36) and the malformed count (h = 4: 440 of 21 045; h = 5: 0) myself. Decision: **accept**. Promoted as **E-301**; **R-015** written; E-263 carries a correction line; H-026 restated (central-only formula falsified, general fit as the new wording, still SUPPORTED at key level).
- **Maverick (T6)** -- purity of key cells equals the shuffle, so "key irrelevant to defect". Referee scholar: **major revision** (purity is blind because both key classes have mode 2; MI p about 0.0005; not-a-function already in E-283; "E-285 candidate" label collides with E-285). Decision: **revise**; onto *Awaiting revision*. The run is recorded as **E-302** (note, not a finding). T6's closure is deferred.

## Consequences
- H-026's central-only formula and the E-263/E-296 "messy (4,7)" clause fell (R-015). E-296's remark "non-lone (4,7)-type classes: messy" should read "(4,7) is a lone-3 class with K0 = 4"; E-296/E-299 are otherwise unchanged. Thread S-1 updated.
- Any other round-097 hit file with entries >= 10 shares the writer bug; none other is cited by the record that I know of.
- The R-candidate "a key-guarded walk stays in one class" does not move: 13 of 13 YES-side chains, still no NO-capable invariant at n <= 7. The F-candidate H-026 blocker "tail ends directly" is met; "a non-lone class" and class level remain.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- Round 113, question 1 (PDFs): citations stay UNVERIFIED -- decided by the chair of round 114; no answer from the human. No overnight run proposed.
