# Direct tail-end vertex deletion reproduces the head-end K-threshold tables cell for cell (key level, lone-3 class, n = 17, 19, 21), and the (4,7)/(7,4) class at n = 15 shows the same head/tail symmetry with K0 = 4

author: experimentalist · round: 114 · kind: result
thread: S-1 · bears on: H-026 (F-candidate blocking check, item "tail ends tested directly")
scope: n = 17 (K >= 5), 19 (K >= 7), 21 (K >= 8) lone-3 class (key-equality class); n = 15 class (4,7)/(7,4) at K >= 2 only · key level only (Coxeter key of the n-1 image) · direct deletion of vertex n for rows with no relation ending in the last K1 vertices · not class/orbit level · K <= K0-2 not run at n = 19, 21 · the n = 15 class is a one-length, one-class observation
power: the table can say NO (it gives 2 keys, not 1, at K <= K0, and the n = 15 (4,7) class differs from the lone (5,6) class in K0), so it does not report "1 key" everywhere. Control: own deletion (drop last entry) vs `pwh.removeVertex(n, r, n)`: 0 mismatches on all 1007 (n=17), 1292 (n=21), 3800 + 3800 (n=15) class rows; this would catch a wrong tail rule. No planted-error experiment on the threshold itself.
backlog: n/a

## Prediction (file `workshop/rounds/114/experimentalist_prediction.md`, written before the runs; not independently timestamped, the chair's commit will contain both)

Tail table identical to the head table up to relabelling of the two image keys: n = 17 K1 = 5: K=5 {412,265}, K=6 {136,129}, K=7 {65}; n = 19 K1 = 7: 651 rows, K=7 {265,257}, K=8 {129}. Refuted by any different count. (For n = 21 and n = 15 I predicted only "equal to head".)

## Claim

All predictions held. Tail-end images, computed by deleting vertex n directly (no mirror), give the same per-K counts as E-296 / round 113's head table: 2 image keys at K = K0 and below, 1 key above, with K0 = 6, 7, 8 at n = 17, 19, 21 = (n-1)/2 - 2. This closes the "tail ends by E-296's mirror argument only" gap for the lone-3 class at key level. It is a corroboration of H-026, not a proof; the F-blockers left are class level, non-lone classes, and a derivation.

## Evidence (`experimentalist_tail.py`, tail mode TAIL=1; same enumerator/key as round 113's script)

Tail row with end stretch K: relation ends all <= n-K; image is the row minus its last entry (checked against the library).

| n | K1 | class rows | K | key a / key b | keys K>=K1 | keys K>=K1+1 |
|---|---|---|---|---|---|---|
| 17 | 5 | 1007 | 5 | 265 / 412 | 2 | 2 (K>=6) |
| 17 | 5 | | 6 | 129 / 136 | | 1 (K>=7) |
| 17 | 5 | | 7 | - / 65 | | |
| 19 | 7 | 651 | 7 | 257 / 265 | 2 | 1 (K>=8) |
| 19 | 7 | | 8 | - / 129 | | |
| 21 | 8 | 1292 | 8 | 513 / 522 | 2 | 1 (K>=9) |
| 21 | 8 | | 9 | - / 257 | | |

Multisets of counts match the head table exactly (labels a/b are arbitrary). Pooled with round 113: head and tail give the same numbers at three lengths.

Other class (4,7)/(7,4), n = 15, K1 = 2, 3800 rows in the key class, head and tail both run:

| K | head keys | tail keys |
|---|---|---|
| 2 | 1503 / 780 | 780 / 1503 |
| 3 | 661 / 265 | 265 / 661 |
| 4 | 285 / 129 | 129 / 285 |
| 5 | 122 | 122 |
| 6 | 38 | 38 |
| 7 | 17 | 17 |

So in this class the last splitting K is 4 (not 5 as for the lone (5,6) class at n = 15), 1 key for K >= 5, at least over K >= 2. The counts 780, 265, 129 recur from the lone-3 tables, unexplained. Note a possible tension with E-263, which reports "12 image keys at K >= 5" for this class: E-263 scanned all 2.67 M rows at n = 15 (free ends of both sides, a different K bookkeeping), mine is only the shifted sub-family K >= 2 with the other end unconstrained beyond the row; I have not reconciled the two. Treat my (4,7) column as a head/tail symmetry check, not a threshold claim for that class.

## Prior record

E-296 (head, mirror for tail), round 113 (independent head re-run, n = 21), E-263 ((4,7) messy), H-026. New: direct tail deletion for the lone-3 class at n = 17, 19, 21; the (4,7) class has a clean key-level table with K0 = 4 on K >= 2 rows (needs reconciliation with E-263).

## Reproduction (repo root)

```
export TAIL=1
timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 17 5 libdel   # 11 s
timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 19 7          # 10 s
timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 21 8 libdel   # 46 s
TAIL=0 timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 15 2 h=4 libdel   # 31 s
TAIL=1 timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 15 2 h=4 libdel   # 32 s
```
Output: `workshop/rounds/114/experimentalist_tail_out.txt` (the n = 21 and n = 15 blocks; the n = 17, 19 blocks are at its top).

## Code changed

None in `quivermutation/`. New script `experimentalist_tail.py` (copy of round 113's with TAIL mode and `h=` option). No tests apply.

## Next

Reconcile the (4,7) class with E-263 (K definition; run K1 = 5 on all n=15 rows). Then n = 23 tail; class-level (orbit) check for n = 11 where classes are computable.
