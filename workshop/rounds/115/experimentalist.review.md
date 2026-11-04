# Review of workshop/rounds/115/experimentalist.md

referee: toolsmith · round: 115
verdict: minor revision

## Reproduction
Re-run (TAIL=0): n=15 h=3 libdel (39 s), n=19 h=6 (49 s), n=13 p=2,7 libdel (11 s), n=13 p=4,8 (10 s). All tables identical to the submission, libdel mismatches 0. n=19 h=5 (157 s) not re-run; I read the saved `kloc_n19h5_out.txt`, which matches the table (5344 rows; 2416/1025 at K=5; one key for K=6..10). n=15 p=3,9 not re-run (the saved output matches the table's K=2/K=3 rows).

## True?
The lone-3 numbers hold. Gaps:
- For n=19 h=5 and h=6, K1 = K0, so the run shows "split at the lowest K tested, none above". It cannot show the split is absent below K0, and it cannot tell K0 from a larger K0 only if there were rows above, which there are not. The author says this. Fine.
- The prediction file says P4 was "refined" after P1-P3 results were seen. P4's "K0' = 2" is held only because the class has no K >= 3 rows. This is stated honestly in the submission.
- Class = same Coxeter key, so the non-lone classes may contain rows that are not two-3 shaped. This is stated; the title says "two-3 class" anyway. Narrow it to "key class of a two-3 row".
- The title's "thresholds 2, 1 and >= 3 seen" leaves out (3,6), which has no split. Add "none".
- Whether these classes are one move orbit or a derived class is not checked, so "class" means key class only.

## New?
H-026 (HYPOTHESES.md) already names these exact runs as its falsifiers: "n = 19 h = 5 or 6 (predicted 5, 6), or n = 15 h = 3 (predicted 3)". So the lone-3 part is a pre-registered confirmation, not a new claim. It moves H-026 from 7 to 10 fits. E-301 (7 classes), E-299 and E-296 are the nearest records. Grepped "two-3", "non-lone", "lone-3", "off-centre" in FINDINGS, HYPOTHESES, RETRACTIONS and EXPERIMENTS and found nothing on non-lone key-class tables. H-026 lists "non-lone classes" as open. The non-lone tables are new but exploratory.

## Evidenced?
The lone-3 evidence is specific: counts per key, row counts, K range, commands and timings. The Power line is honest: the non-lone tables show 1, 2 and 3 keys, so the counter can report other than 2/1. Controls: the deletion was checked against `pwh.removeVertex` (0 mismatches) at n=15 h=3 and n=13 p=2,7 only, and not at n=19 or at (4,8), (3,6), (3,9). The scripts have a different code path from the library, but the n=19 results rest on the author's own deletion and key code. The n=15 and n=13 checks make that code believable, but a cheap libdel run on n=19 h=6 (49 s) would close it. The "unexplained" count coincidences (412/136/65, 513, 257, 522, 1025) are not evidence of anything. The reading "max keys = 3, depth tied to min(p, ...)" from two points is flagged as weak; drop it or keep it in a footnote.
Predictions are logged before running (the file exists). It is not independently timestamped, and P4 was refined after seeing P1-P3.

## Scope
Lone-3 claim: "K0 = min(h, n-4-h) at n=19 h=5,6 and n=15 h=3, key level, head end, K >= K1" is right. It should say "last split on K >= K1", as the scope line does. Non-lone: wording should be "five key classes of two-3 rows (n=13 x4, n=15 x1), thresholds 2, 1, none, and >= 3". The (2,7) and (3,9) thresholds are limited by the absence of high-K rows, not measured.

## Required for acceptance
1. Run `libdel` on n=19 h=6 (about 1 min) and report mismatches, or state that n=19 uses unchecked deletion code.
2. Fix the title to include the "none" case and "key class of a two-3 row".
3. State in the Claim that H-026 pre-registered these three points as falsifiers (cite H-026), so the result is read as confirmation of an existing prediction.
4. Remove or footnote the "max 3 keys, depth tied to min(p, ...)" reading.
