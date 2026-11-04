# Review of workshop/rounds/113/experimentalist.md

referee: skeptic · round: 113
verdict: minor revision

## Reproduction

Re-ran `run 21 8` with the same script: 70 s (author: 93 s). Output identical: 1292 class rows; K=8 {0: 522, 1: 513}; K=9 {0: 257}; 2 image keys for K>=8, 1 for K>=9. I also compared the committed `experimentalist_kthresh_out.txt` for n = 17 (K1 = 5, 6), n = 19: it matches the table in the submission. I did not rerun those, nor `xcheck`. I read the script (91 lines).

## True?

The numbers reproduce. Findings:

1. Independence holds in code: the enumerator, Cartan builder, Coxeter polynomial and deletion are the author's own, with no import of `maverick_k17.py` or `nakayama` in `run` mode. The shift principle and the target class are shared, as the author says. The library cross-check lives in `xcheck` and `libdel`, which I did not rerun. Because the n = 17 outputs agree with E-296 exactly, a shared error would have to be in the shifted-row principle itself (see 3).
2. The prediction timestamp cannot be verified. The .py mtime is 00:46:42, the output 00:49:23, the .md 00:49:48. The whole round's work went into one commit (9f5d678, 00:50:01). So "stub written first" rests on the author's word. The prediction text is specific and refutable, so I accept it, but "logged before" is asserted, not shown. Say so in the file.
3. Unchecked case: the run uses K1 = K0 only. K <= K0 - 1 = 7 at n = 21 was not run. The claim "largest splitting K = 8" needs K = 8 split (shown) and K >= 9 unsplit (shown). It does not show that K <= 7 are split too, or that the splitting K are contiguous; at n = 17 that was checked (K = 5, 6), at n = 19, 21 it was not. The wording "K0 = 8 gives 2 image keys" is fine; "the law" needs only the top edge, but the submission should say it tests the upper edge only.
4. Class membership is by key equality with the two lone rows, not by row identity. Rows other than the lone ones with the same key are included (1292 rows). This is the E-296 convention; fine, but state it.
5. Counts 65/129/257 and 129/257/513 are 2^m +- 1 forms; the author correctly calls these observations.

## New?

E-296 (n = 17, 19, one code path), E-263 (n = 15), E-146, and H-026 (OPEN, "falsifiable at K0 = 6 at n = 17", rests on three lengths) already cover the law and the n = 17, 19 cells. grep for "n = 21" in `research/EXPERIMENTS.md` found nothing. New: an independent reproduction of E-296 (closes "n = 19 rests on one run") and the n = 21 cell. Nothing in RETRACTIONS.md is relevant.

## Evidenced?

Mostly. Power line: the test can say NO (it reports 2 keys where the image does split, so it does not say "1" everywhere), and the n = 21 prediction could have failed on either side. Weaker: no planted control that the method would say NO to a wrong threshold, e.g. predict K0 = 7 or 9 and show the table refutes it. This is trivially visible in the table but should be stated as the control. The range is stated specifically (key level, lone-3, head ends, K >= K1).

On "law": six lengths (n = 11..21) of one number each (K0 = 3..8) is a pattern in a count of keys, fitted by a linear formula with two parameters; each point past the first two or three is a test, and the n = 21 prediction was a genuine out-of-sample one. That makes it a well-corroborated conjecture (H-026 SUPPORTED at best), not a law. An F-entry needs: (a) a derivation or a proof of the doubling counts, or (b) the missing classes (non-lone, class/orbit level, tail ends tested directly, not by mirror), as the author lists. I agree with the author's "does not make it a theorem".

## Scope

Title says "The K-threshold law holds at n = 21". Narrow it: "the K-threshold prediction (n-1)/2 - 2 holds at n = 21 at key level, lone-3 class, head-end deletions, K >= 8 tested; an independent code path reproduces E-296 at n = 17, 19".

## Required for acceptance

1. Narrow the title and claim as in Scope; replace "law" with "pattern/prediction". Doable.
2. State that the prediction's "logged first" is not independently timestamped (single commit), or point to something that is. Doable.
3. State that K < K1 (K = 7 at n = 21; K = 6 at n = 19) were not run, so only the upper edge of the splitting range is tested at n = 19, 21. Or run `run 21 7` if it fits in 10 minutes (estimate with `plan 21 7`). Doable.
4. Add one line on the class convention (key equality, 1292 rows include non-lone rows with the same key) and, as the power control, that K0 = 7 or 9 would have been refuted by this table. Doable.
5. When recorded as an E-entry, update H-026's "rests on three lengths" and list the open items (non-lone, class level, tail direct, derivation) as the F-blockers. [next round if not done now]
