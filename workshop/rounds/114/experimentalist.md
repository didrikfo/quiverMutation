# Direct tail-end vertex deletion reproduces the head-end K-threshold tables cell for cell (key level, lone-3 class, n = 17, 19, 21); the K-threshold of a lone-3 key class is min(h, n-4-h), so H-026's (n-1)/2 - 2 holds for the central class only; E-263's (4,7) line is a data artefact

author: experimentalist · round: 114 (revised after referee) · kind: result
thread: S-1 · bears on: H-026 (F-candidate blocking check, item "tail ends tested directly"; falsifier "a lone-3 row whose threshold differs")
scope: central lone-3 class at n = 17 (K >= 5), 19 (K >= 7), 21 (K >= 8) · non-central lone-3 classes (h < central) at n = 13 (h = 0..4, K >= 2, head and tail) and n = 15 (h = 4, K >= 2), n = 17 (h = 5, K >= 3), tail rows only for K above the stated K1; "K0 = h" is read off K >= K1 tables, K <= K1-1 not run · key level only (Coxeter key of the n-1 image) · direct deletion of vertex n · not class/orbit level · K <= K0-2 not run at n = 19, 21 · n = 19, 21 rows were NOT re-run by the referee (author's run only; the n = 21 table has only K = 8, 9)
power: the tables can say NO (2 keys at K <= K0, and non-central classes give a different K0 from the central one, so the table does not report "1 key" everywhere). Control: own deletion (drop last entry) vs `pwh.removeVertex(n, r, n)`: 0 mismatches on every class row of every run listed (1007, 1292, 3800 + 3800, 6600, and the n = 13 runs). No planted-error experiment on the threshold itself.
backlog: n/a

## Response to referee (theorist, minor revision)

1. **E-263 (4,7)/(7,4) "12 image keys at K >= 5" explained.** Cause found, and it differs slightly from the review's. Unpack and count (repo root):
   ```
   mkdir -p /tmp/h && tar xzf workshop/rounds/097/experimentalist_scan15_hits.tar.gz -C /tmp/h
   .venv/bin/python workshop/rounds/114/experimentalist_hitscount.py /tmp/h     # output: experimentalist_hitscount_out.txt
   ```
   h = 4 files: 21 045 lines (= E-263's class size), 20 605 of length 13, 374 of length 14, 66 of length 15 (440 malformed). The malformed lines are NOT shard-boundary joins (the review's guess; `cat` of the files does join them, per-file reading does not): `experimentalist_scan15.py` line 36 writes `''.join(map(str, r))`, so an entry >= 10 is two characters, and `ends` (line 43) parses `int(c) for c in l` digit by digit, giving wrong rows with K = -2, -1 or large. Check (`experimentalist_hitsvs.py`, output `_hitsvs_out.txt`): all 3800 tail-K>=2 class rows of my enumerator appear in the files as strings, 44 of them with an entry >= 10 (14 or 15 characters); so no class row is "lost" (the review's 3756 vs 3800 and "44 absent" is the same 44, present but unparsable). Correct reading: with a separator-safe parse, K=5 {122}, K=6 {38}, K=7 {17}, one key each (the review's table), K0 = 4. **Proposed E-263 correction:** the (4,7)/(7,4) clause "messier (12 image keys at K >= 5)" is wrong; replace by "K0 = 4: K = 3, 4 give 2 image keys, K >= 5 give 1" (E-263's own K >= 3 window), and add a RETRACTIONS entry for that clause if the chair agrees; the (5,6) clause stands (item 3). E-296's "non-lone (4,7)-type classes (E-263: messy)" should read "(4,7) is a lone-3 class with K0 = 4". The same writer bug (digit join) is a warning for any other round-097 hit file with entries >= 10.
2. **H-026 falsifier and the formula.** New runs (`experimentalist_tail.py`, `h=`; outputs `_n13_out.txt`, `_n17h5_out.txt`, `_n17h5head_out.txt`): last K with 2 image keys, on K >= K1 rows, head = tail up to key labels:

   | n | class (h, n-4-h) | rows | K0 observed | central (n-1)/2-2 |
   |---|---|---|---|---|
   | 13 | (4,5) | 788 | 4 | 4 |
   | 13 | (3,6) | 880 | 3 | 4 |
   | 13 | (2,7) | 1150 | 2 | 4 |
   | 13 | (1,8), (0,9) | 2087, 576 | no split at K >= 2 (K0 <= 1) | 4 |
   | 15 | (5,6) | 18 328 (E-263) | 5 | 5 |
   | 15 | (4,7) | 3800 | 4 | 5 |
   | 17 | (6,7) | 1007 | 6 | 6 |
   | 17 | (5,8) | 6600 | 5 (K = 3, 4, 5 split, K >= 6 one key) | 6 |

   So H-026 as worded ("a lone-3 row at n = 11..15 whose threshold differs") **is falsified**: the n = 13 classes h = 3, 2 and n = 15 h = 4 have thresholds 3, 2, 4. The (n-1)/2 - 2 formula survives only for the central class (h = (n-4)/2 rounded down, i.e. both lone-3 rows nearest the middle). On all 7 classes with a split the data fit K0 = min(h, n-4-h) = the shorter free stretch beside the lone 3 (n = 19, 21 central: 7, 8, also consistent). This is a fit to seven non-trivial points, not a theorem, and K0 <= K1-1 classes (h = 0, 1 at n = 13) say only "no split at K >= 2". Suggested H-026 rewrite: K0 = shorter free stretch, central class is the special case; the falsifier becomes "a lone-3 class with K0 != min(h, n-4-h)".
3. **E-263 (5,6) class.** h = 5 files: 18 328 lines, all of length 13, 0 malformed (`_hitscount_out.txt`); so no row of the (5,6) class is damaged, and no malformed row can sit among its K >= 3 ends. E-263's (5,6) numbers are unaffected by this defect. (E-296's agreement at K = 5, 6 was consistent with this.)
4. Dropped: "not a threshold claim", "not reconciled", "needs reconciliation". The (4,7) K0 = 4 claim is supported on tail K >= 2 rows (the range stated in scope); the head columns at n = 15 are a second run of the same enumerator, not independent of E-263's data (but see item 1: E-263's data are corrupted for this class).
5. The n = 19 and n = 21 rows and the TAIL=0 n = 15 run were not re-run by the referee (it read the blocks in `_tail_out.txt`); the n = 21 table covers K = 8, 9 only, so "cell for cell" is for the cells E-299 has. Prediction file covered n = 17, 19 only.

## Claim

All round-114 predictions held. Tail-end images computed by deleting vertex n directly give the same per-K counts as E-296's head table: 2 image keys at K = K0 and below, 1 key above, K0 = 6, 7, 8 at n = 17, 19, 21 (central class). New: the lone-3 key classes off the centre have smaller K0 (min(h, n-4-h)), which falsifies H-026's wording and restricts its formula to the central class; and E-263's (4,7) "messy" line is a digit-join writer artefact. A corroboration with a correction, not a proof; left: class level, non-lone classes, a derivation of K0 = shorter stretch.

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

So in this class the last splitting K is 4, 1 key for K >= 5, over K >= 2 (E-263's 12-key line is the artefact of response item 1). The counts 780, 265, 129 recur from the lone-3 tables, unexplained.

## Prior record

E-296 (head, mirror for tail), round 113 (independent head re-run, n = 21), E-263 ((4,7) "messy": wrong, response item 1), H-026. New: direct tail deletion at n = 17, 19, 21; clean K0 tables for the non-central lone-3 classes (response item 2).

## Reproduction (repo root)

```
export TAIL=1
timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 17 5 libdel   # 11 s
timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 19 7          # 10 s
timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 21 8 libdel   # 46 s
TAIL=0 timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 15 2 h=4 libdel   # 31 s
TAIL=1 timeout 10m .venv/bin/python workshop/rounds/114/experimentalist_tail.py run 15 2 h=4 libdel   # 32 s
```
Reply runs: `TAIL=0/1 ... run 13 2 h=H libdel` (H = 0..4, 3 s each), `TAIL=1/0 ... run 17 3 h=5 [libdel]` (123 s).
Output: `workshop/rounds/114/experimentalist_tail_out.txt` (the n = 21 and n = 15 blocks; the n = 17, 19 blocks are at its top).

## Code changed

None in `quivermutation/`. New scripts `experimentalist_tail.py` (copy of round 113's with TAIL mode and `h=` option), `experimentalist_hitscount.py`, `experimentalist_hitsvs.py`. No tests apply.

## Next

Chair: decide the E-263 correction / RETRACTIONS entry and the H-026 rewrite. Then: test K0 = min(h, n-4-h) at n = 19 (h = 5, 6) and n = 15 h = 3 (K1 = 2); n = 23 tail; class-level check at n = 11.
