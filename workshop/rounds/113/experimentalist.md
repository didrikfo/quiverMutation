# The K-threshold law holds at n = 21 (K = 8 gives 2 image keys, K >= 9 gives 1) and the E-296 results at n = 17, 19 reproduce exactly by an independent code path

author: experimentalist · round: 113 · kind: result
thread: S-1 · bears on: H-026 (F-candidate blocking check)
scope: n = 17, 19, 21 · lone-3 key class only ((6,7)/(7,6), (7,8)/(8,7), (8,9)/(9,8)) · key level only (Coxeter polynomial of the n-1 image) · head-end deletions only (tail by E-296's mirror argument, not re-tested) · no depth/time cap, all rows with a head end >= K1 enumerated (K1 = 5/6 at n = 17, 7 at n = 19, 8 at n = 21; n = 21 K1 = 9 not run, it is a subset of K1 = 8) · not class level, not non-lone classes
power: the test does say NO (image keys > 1 = "not transported") at K <= K0 on the E-263 / E-296 cells where that is known: it reports 2 keys at n = 17 K = 5, 6 and n = 19 K = 7, i.e. it does not report 1 key everywhere. No planted-error control; the cross-check that would have caught a code error is the library comparison below (0 mismatches).
backlog: n/a

## Prediction, logged in this file before any n = 21 run

Law K0 = (n-1)/2 - 2 gives K0 = 8 at n = 21. Predicted: lone-3 class (8,9)/(9,8); head-end images at n = 20 have 2 keys for K >= 8 and exactly 1 key for K >= 9, with K = 8 the only splitting stretch. Refuted by: 1 key at K >= 8 (threshold lower, K0 <= 7), or 2 keys at K >= 9 (K0 >= 9). (The stub was written to this file first; the command history shows the n = 21 run after it.)

## Claim

Independent re-run agrees with E-296 on every cell, and the n = 21 prediction was confirmed. Largest splitting K = 6 (n = 17), 7 (n = 19), 8 (n = 21), the law (n-1)/2 - 2, now at a sixth length (3, 4, 5, 6, 7, 8 at n = 11..21). Key level, lone-3 class, head ends only; says nothing about class level, orbits, non-lone classes, or tail ends beyond E-296's mirror test. A disagreement would have looked like: a different count of rows in the key class (E-296: 330 at n = 17 K1 = 6), a different image-key number at a K cell, or the 265/129/65 pattern failing.

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

Blocking checks remaining for an F-status of H-026: (1) the non-lone classes (E-263's (4,7)-type is messy), (2) class level/orbit level rather than key level, (3) tail ends beyond the mirror test (the same script with a tail shift would be an independent re-test), (4) a derivation of 2^m + 1 / 2^m - 1 counts. Theorist: why the law is (n-1)/2 - 2. Experimentalist next: n = 23 K1 = 9 (est. ~6 min; shifted length 14 = 742 900 rows, plan first).
