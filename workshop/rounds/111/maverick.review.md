# Review of workshop/rounds/111/maverick.md

referee: theorist · round: 111
verdict: minor revision

## Reproduction

Re-run: `check` (seconds): head sets equal (1430, 429), tail not a shift (13260, 6188), as stated. `run 15 4` (66 s): K4 136/222, K5 65/71, K6 0/33, None 0, image keys K>=5: 2, K>=6: 1, matching the table and half of E-263 (142/130, 66/0, labels swapped). `run 17 6` (94 s): 330 in class, K6 129/136, K7 0/65, keys K>=6: 2, K>=7: 1, matching. Not re-run: `run 17 5` (370 s), `run 19 7` (518 s, near the limit); n = 19 is therefore unverified by me, and it is the claim in the title.

## True?

The shift reduction is sound (head sets equal, shown at n = 12 only, K1 = 3, 4; it is an obvious bijection, so fine). The mirror-closure assumption, which I tried to break: on 400 random n = 11 rows, `lnaCoxeterKey` is mirror-invariant (0 failures), and `removeVertex(n, r, n)` vs `removeVertex(n, mirror(r), 1)` give the same image key (0 failures, 0 None). So tail ends add no new image keys, given the sample. The assumption is plausible but is only argued in the submission from E-263 totals at n = 15; it is not tested at n = 17/19, and E-173 records head and tail images of one row differing at n = 17, which the text does not reconcile (differ as rows, equal as keys is what the above suggests). One class-level slip: table column "in key class" for n = 15 says 358, but the run gives 527 (358 is the K = 4 row count; 136 + 33 the rest); n = 17 "1 007" and n = 19 "651" likely have the same K1-versus-total mix-up (n = 17 K1 = 6 gives 330 = 129+136+65, consistent with total, so 1 007 vs 265+412+... = 412+... should be checked). Wording "K0 = 6 at n = 17" is used in the opposite sense from H-026 ("K >= K0 transports"; H-026 K0 = (n-1)/2-2 = 6 and E-263 K0 = 5 at n = 15 where K = 5 splits). The paper's own gloss flips between the two readings; with the author's convention the statement "K >= 7 transports at n = 17" is what was verified. Fix the wording.

## New?

Nothing found in FINDINGS/RETRACTIONS for the n = 17/19 K table. Prior: E-146, E-173, E-263, H-026 (its falsifier "K0 = 6 at n = 17" is what is tested; not falsified). New: n = 17 and n = 19 whole-class key-level table, and the shift reduction (cheap method; E-263's note on two-pass plan is the nearest).

## Evidenced?

Mostly. Counts per cell are given, control reproduces known data (so the method does say "two keys" where known). Power line is honest: no planted error; but the control tests the same code that could share a bug. Weak points: (a) n = 19 covers K >= 7 only, so it shows K = 7 splits and K = 8 transports, but "K >= 7 doesn't transport" does not mean K = 6 was seen; the lower side (K0-1 splits) is the K = 7 cell, fine. (b) "K = 8 : 129/0": one image key, but this is a single class cell, no check the K = 8 key equals the K >= 7 key at n = 18 across other K (only 1 value of K). (c) the class is the key class of the lone-3 rows h = 7, 8 only, and the (4,7)-type messy classes (E-263) are not examined, so a law of "the lone-3 key class" is claimed for one class per n. (d) the regularity 265/129/65 is observed, unexplained, flagged honestly.

## Scope

Title says "K0 = 6 holds at n = 17 and K0 = 7 at n = 19" with qualifiers key level, one class, head ends. Narrowed wording: "for the single key class containing the middle lone-3 rows, head-end deletions (tail by mirror, sampled equivariance at n = 11 only), K >= K0+1 gives one n-1 key and K = K0 gives two, at n = 17 (verified by me) and n = 19 (not re-run)". "Fourth and fifth length" is fair for key level.

## Required for acceptance

1. Fix the "in key class" column (527 / n = 17 K1=5 / n = 19 totals) and state which count each is.
2. State the K0 convention once and align with H-026 ("K >= K0 transports" vs "largest splitting K").
3. Add the mirror test (key mirror-invariance and removeVertex equivariance, my script: random n = 11 rows, 400, 0 failures) or one n = 15 full-class check, as the basis for dropping tail ends; one sitting.
4. Reconcile with E-173 (head vs tail images of one row differ as rows; same key).
5. Narrow the title or say "one key class per n".
