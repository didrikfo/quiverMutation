# The K-threshold prediction (n-1)/2 - 2 holds at n = 21 (K = 7, 8 give 2 image keys, K >= 9 gives 1; key level, lone-3 class, head-end deletions) and E-296 at n = 17, 19 reproduces by an independent code path

## Response to referee

Items 1-4 done; item 5 is for next round (it applies when this is recorded as an E-entry: update H-026's "rests on three lengths", list the F-blockers).
1. Title/claim narrowed; "law" replaced by "pattern/prediction" throughout (Claim, Prior record, Next). Status: a corroborated prediction, H-026 SUPPORTED at best, not a law.
2. Agreed. The "logged first" prediction is NOT independently timestamped: the whole round is one commit (9f5d678), and file mtimes (.py 00:46, out 00:49, .md 00:49) do not show order. It rests on the author's word; the text is specific and refutable, nothing more is claimed.
3. Lower-edge cells run (plan 21 7 est. 264 s, plan 19 6 est. 47 s). Actual: `run 19 6` 54 s, 1953 class rows, K=6 {0: 780, 1: 522}, K=7 {0: 265, 1: 257}, K=8 {0: 129}; 2 keys for K>=6 and K>=7, 1 for K>=8. `run 21 7` 253 s, 3829 class rows, K=7 {0: 1502, 1: 1035}, K=8 {0: 522, 1: 513}, K=9 {0: 257}; 2 keys for K>=7 and K>=8, 1 for K>=9. So K = 7 (n = 21) and K = 6 (n = 19) split too; with K1 = 8 / 7 these are the cells just below K0, so splitting is contiguous over the tested K = K0-1, K0 (n = 21: 7, 8; n = 19: 6, 7), then 1 key above. Not run: K <= K0-2 (n = 21: 6 and lower; n = 19: 5 and lower); contiguity below that is unchecked. The K=8 cell is the same 522/513 as in the K1 = 8 run (shift principle is consistent across K1).
4. Class convention: class membership is by key equality with the lone rows' key, not row identity; the 1292 rows (K1 = 8) include non-lone rows with the same key (the E-296 convention). Power control: the table would have refuted a wrong threshold. K0 = 7 predicts 1 key at K >= 8: the table shows 2 (522/513), refuted. K0 = 9 predicts 2 keys at K >= 9: the table shows 1 (257), refuted. Likewise at n = 19 (K0 = 6 refuted by 2 keys at K = 7; K0 = 8 by 1 key at K = 8). This is a check by reading the table, not a separately run planted-error experiment.


author: experimentalist · round: 113 · kind: result
thread: S-1 · bears on: H-026 (F-candidate blocking check)
scope: n = 17, 19, 21 · tested K range: n = 21 K = 7, 8, 9+; n = 19 K = 6, 7, 8+; n = 17 K = 5, 6, 7+ (K <= K0-2 at n = 19, 21 not run) · lone-3 key class only ((6,7)/(7,6), (7,8)/(8,7), (8,9)/(9,8)) · key level only (Coxeter polynomial of the n-1 image) · head-end deletions only (tail by E-296's mirror argument, not re-tested) · no depth/time cap, all rows with a head end >= K1 enumerated (K1 = 5/6 at n = 17, 6/7 at n = 19, 7/8 at n = 21) · not class level, not non-lone classes
power: the test does say NO (image keys > 1 = "not transported") at K <= K0 on the E-263 / E-296 cells where that is known: it reports 2 keys at n = 17 K = 5, 6 and n = 19 K = 7, i.e. it does not report 1 key everywhere. Power control by table reading: K0 = 7 or 9 (n = 21) would have been refuted (see Response 4); no separately run planted-error experiment; the cross-check that would have caught a code error is the library comparison below (0 mismatches).
backlog: n/a

## Prediction, written in this file before the n = 21 run (not independently timestamped, single commit)

Pattern K0 = (n-1)/2 - 2 gives K0 = 8 at n = 21. Predicted: lone-3 class (8,9)/(9,8); head-end images at n = 20 have 2 keys for K >= 8 and exactly 1 key for K >= 9, with K = 8 the only splitting stretch. Refuted by: 1 key at K >= 8 (threshold lower, K0 <= 7), or 2 keys at K >= 9 (K0 >= 9). (The stub was written to this file first per the author; this is unverifiable, see Response 2.)

## Claim

Independent re-run agrees with E-296 on every cell, and the n = 21 prediction was confirmed, including K = 7 at n = 21 and K = 6 at n = 19 (both 2 keys). Largest splitting K = 6 (n = 17), 7 (n = 19), 8 (n = 21), the pattern (n-1)/2 - 2, now at a sixth length (3, 4, 5, 6, 7, 8 at n = 11..21). Key level, lone-3 class, head ends only; says nothing about class level, orbits, non-lone classes, or tail ends beyond E-296's mirror test. A disagreement would have looked like: a different count of rows in the key class (E-296: 330 at n = 17 K1 = 6), a different image-key number at a K cell, or the 265/129/65 pattern failing.

## Evidence

Different code path: own row enumerator (DFS on the antichain condition; count equals Catalan numbers, equals `nakayama.allRelationLengths` for m = 5..11), own Cartan builder, own Coxeter polynomial (numpy charpoly of -C^{-T}C, rounded, rounding residual asserted < 1e-6), own vertex-1 deletion (drop the first entry). Shared with maverick: only the shift principle (row with head end >= K1 = 0^K1 + arbitrary row) and the choice of target class. Cross-check against the library: key 0 mismatches on 1200 random rows n = 9..11; deletion (drop first entry) vs `pwh.removeVertex(n, r, 1)` 0 mismatches on 1859 rows (n = 9, 10) and on all 330 class rows at n = 17 K1 = 6.

Image keys at n - 1, rows by end stretch K (two image keys: counts per key):

| n | K1 run | class rows | K | key 0 / key 1 | image keys K>=K1 | K>=K1+1... |
|---|---|---|---|---|---|---|
| 17 | 6 | 330 | 6 | 136 / 129 | 2 | K>=7: 1 |
| 17 | 5 | 1007 | 5 | 412 / 265 | 2 | K>=6: 2, K>=7: 1 |
| 17 | (5) | | 6 | 136 / 129 | | |
| 17 | (5) | | 7 | 65 / - | | |
| 19 | 7 | 651 | 7 | 265 / 257 | 2 | K>=8: 1 |
| 19 | 7 | | 8 | 129 / - | | |
| 21 | 8 | 1292 | 8 | 522 / 513 | 2 | K>=9: 1 |
| 21 | 8 | | 9 | 257 / - | | |

Agreement with E-296: n = 17 K = 5 two keys, K = 6 two, K >= 7 one; n = 19 K = 7 two, K = 8 one: exact. The single-key cell sizes 65, 129, 257 and the per-K counts reproduce the doubling pattern (65 at n = 17 K = 7, 129 at n = 19 K = 8, 257 at n = 21 K = 9; split cell's key-1 count 129, 257 513 = 2^m + 1 style) noted as unexplained in E-296; it is an observation, not a derivation. The n = 19 K = 7 cell, which E-296 said rested on one run, is now confirmed by a second implementation.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/113/experimentalist_kthresh.py xcheck          # ~10 s
timeout 10m .venv/bin/python workshop/rounds/113/experimentalist_kthresh.py plan 21 8       # 1 s (est. 97 s)
timeout 10m .venv/bin/python workshop/rounds/113/experimentalist_kthresh.py run 17 6 libdel # 8 s
timeout 10m .venv/bin/python workshop/rounds/113/experimentalist_kthresh.py run 17 5         # 21 s
timeout 10m .venv/bin/python workshop/rounds/113/experimentalist_kthresh.py run 19 7         # 24 s
timeout 10m .venv/bin/python workshop/rounds/113/experimentalist_kthresh.py run 21 8         # 93 s
```
Outputs: `workshop/rounds/113/experimentalist_kthresh_out.txt` (small).

## Prior record

E-296 (n = 17, 19), E-263 (n = 15), E-146 (n = 11..13), H-026 (falsifiable at K0 = 6 at n = 17; its n = 19 and now n = 21 predictions held). New here: an independent implementation reproduces E-296 exactly (closing "n = 19 rests on one run"), and the n = 21 cell. The n = 21 result extends the pattern to a sixth length at key level; it does not make it a theorem.

## Code changed

None in `quivermutation/`. New script `workshop/rounds/113/experimentalist_kthresh.py`; no tests apply (no library file touched).

## Next

Blocking checks remaining for an F-status of H-026: (1) the non-lone classes (E-263's (4,7)-type is messy), (2) class level/orbit level rather than key level, (3) tail ends beyond the mirror test (the same script with a tail shift would be an independent re-test), (4) a derivation of 2^m + 1 / 2^m - 1 counts. Theorist: why the pattern is (n-1)/2 - 2. Experimentalist next: n = 23 K1 = 9 (est. ~6 min; shifted length 14 = 742 900 rows, plan first).
