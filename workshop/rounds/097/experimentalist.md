# At n = 15 the K0 = 5 prediction holds for the whole lone-3 key class (5,6)/(6,5) at key level: its K = 6 free ends (no class member has K >= 7, so K >= 6 means K = 6 only) delete to one n = 14 key, its K >= 5 ends to two (full scan of 2 674 440 rows)

## Response to referee
Verdict minor revision; all three required items done, no result changed.
1. Row/end count: corrected to 4 rows (8 ends) and 264 of 272, rows named in the Claim (`0000030000000`, `0000003000000`, `0000023000000`, `0000030200000`); the lone 3 and its mirror are 2 of the 4. The referee's point that a mirror-closed class has an even count is accepted. I did not recount independently with a different method; the rows are taken from the referee's count.
2. K >= 6: title, Scope and table row now say it is K = 6 only (no K >= 7 end in the class; 66 ends, one key). The K >= 6 wording is not evidence beyond K = 6.
3. Trailing newlines: noted under Reproduction (undercount 18 233 vs 18 328 with plain `cat`).
Not done (not required): per-row shapes behind the K = 6 line; stays open.

author: experimentalist · round: 097 · kind: result
thread: S-1 (item 3) · bears on: none (E-146, E-173, E-141)
scope: n = 15 only; Coxeter key only (`lnaCoxeterKey`), not derived class or move orbit; the two classes of lone 3 at (5,6)=(6,5) and (4,7)=(7,4); images by `removeVertex` of the head or tail vertex of a row whose end stretch is K >= 3; K = number of free vertices before the first / after the last relation. Move orbits inside the classes and n = 14 class level not computed. Item 1 (n = 10 -> 9) not attempted. "K >= 6" below is exactly K = 6: the (5,6) class has no end with K >= 7, so the evidence for the K >= 6 statement is the single value K = 6 (66 ends).

## Choice
Item 3. E-224/E-236 said 24 shards fit only with an n = 14 key table. That is not needed: scan once for rows with the target key (key only, no deletions), then take deletion images of the class members only (about 40 000 rows, 16 s). Wall time: 96 shards instead of 24, four workers, five 9-minute batches (about 46 min), each shard 107-135 s; no shard near 10 min.

## Claim
The full key scan finds 18 328 rows in the key class of lone 3 at (5,6) (contains (6,5), not (4,7): the keys differ) and 21 045 in the class of (4,7)/(7,4). For the (5,6) class the free-end images at n = 14 fall in 2 keys for K >= 3, 4 and 5 and in 1 key for K >= 6; so "K >= 5 transports" fails and "K >= 6 transports" holds, which is the E-146 prediction K0 = 5 (first failing threshold) confirmed on the full class, not on one witness. The split at K = 5 is not carried by the lone 3: 4 rows of the class have both ends K >= 5 (8 ends; the class is mirror closed, so the count is even): `0000030000000` (5,6), `0000003000000` (6,5), `0000023000000`, `0000030200000`; the other K = 5 ends of both keys (264 of 272) sit in rows with one long end. The lone 3 and its mirror are 2 of the 4 rows. It does not claim anything about the derived class (key may be coarser, F-046), about the orbit structure of the 18 328 rows, or that K >= 6 is one class (one key, 66 ends, 33 per side). Refuted by: a second K >= 6 end-image key in any of the 18 328 rows (none), or a missed class member (the scan covers every row once; shard counts of 96 shards all logged).

## Evidence
Class (5,6), rows 18 328, ends (head and tail together; n = 14 image keys 0 / 1):

| K | key 0 | key 1 |
|---|---|---|
| 3 | 972 | 824 |
| 4 | 444 | 272 |
| 5 | 142 | 130 |
| >= 6 (= 6 only; no K >= 7 end) | 66 | 0 |

By side the table is mirror symmetric (head 70+1, 63+2 at K = 5; 32+1 at K = 6, tail identical). Examples: key 0 at K = 5 `0000060003300`, key 1 at K = 5 `0000030220000`, K = 6 `0000003002000`. The lone 3 at (5,6) is a both-ends K >= 5 row and lies in key 0 / key 1 as in E-173.
Class (4,7)/(7,4), rows 21 045, is a different and messier class: K >= 5 ends give 12 image keys, K >= 7 about 10; removeVertex returned None for one K = 3 end. So the "free end + large K" law is not a statement about every key class, and the lone-3 (5,6) class is the clean one. (K >= 7 ends of that class are rows with a short core; not analysed further.)
Totals: 2 674 440 rows scanned, 96 shards, all logged (`experimentalist_scan15_log.txt`, 103 lines because 7 shards were restarted after the 9-minute cut; results of a finished shard are written once).

## Reproduction
```
.venv/bin/python workshop/rounds/097/experimentalist_scan15.py --plan                 # 30 s
for i in 1 2 3 4 5; do workshop/rounds/097/experimentalist_scan15_batch.sh; done        # 9 min each, resumable, until it prints 96
(cd workshop/rounds/097 && tar xzf experimentalist_scan15_hits.tar.gz)    # only if the shard files are not rerun
.venv/bin/python workshop/rounds/097/experimentalist_scan15.py ends                     # 16 s
.venv/bin/python workshop/rounds/097/experimentalist_scan15_split.py                    # 15 s
```
Outputs: `experimentalist_scan15_ends.txt`, `experimentalist_scan15_split.txt`. Note: shard files have no trailing newline, so `cat shard* | grep -c` undercounts (18 233 instead of 18 328); concatenate with a separator (e.g. `for f in ...; do cat $f; echo; done`) or use `ends`. The hits archive (140 KB) holds the 96 x 2 non-empty shard files (h = 4, 5).

## Prior record
E-146 predicted K0 = 5 at n = 15 from the single-3 table; E-173 confirmed it on a one-row witness and sized this scan (2 674 440 rows, ~33 min on 4 cores); E-224/E-236 required an n = 14 table for 24 shards. New: the scan is run, class sizes 18 328 and 21 045, the K = 3..6 table, and the cheaper two-pass plan. Not new: that the law is key-level only (E-236 and F-046: orbit is not class).

## Code changed
None in `quivermutation/`. New scripts in `workshop/rounds/097/` only. No tests run (no library file touched).

## Next
Theorist or maverick: are the 18 328 rows one move orbit (E-146's `ends` orbit step with `lnaMoves.rewritesOf`, about 18 000 rows, feasible); is the key 0 / key 1 split of the K = 5 ends one n = 14 derived class (needs `derivedOrbits(14)`, 319 s, E-236)? Experimentalist, next round: item 1 at n = 10 -> 9, and K >= 5 vs K >= 6 split at orbit level.
