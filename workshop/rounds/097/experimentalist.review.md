# Review of workshop/rounds/097/experimentalist.md

referee: theorist · round: 097
verdict: minor revision

## Reproduction
- Extracted the hits archive, re-ran `scan15.py ends` (18 s) and `scan15_split.py`: both outputs identical to the committed `_ends.txt` / `_split.txt`. Every K x key number in the table reproduces from the committed shard files.
- Re-ran shard 17 of 96 from scratch (104 s): `_hits_5_17` and `_hits_4_17` byte-identical to the archived ones. The full 46 min scan was not re-run. Log: 103 lines, 96 distinct shards, each "scanned 2674440" (the full enumeration count, so every shard saw every row).
- Class sizes: 18 328 distinct rows (no duplicates across shards) for (5,6); the (4,7) union is the 21 045 stated. A plain `cat | grep -c` gives 18 233 because shard files have no trailing newline. That is an artefact, not an error, but a reader checking by `cat` will be misled.

## True?
The headline reproduces: K >= 6 ends give one n = 14 key, K >= 5 give two (142/130), at key level, for this class.
One stated number is wrong. "Only 3 rows of the class have both ends K >= 5 (6 ends)": I count 4 rows with both ends >= 5: `0000030000000` (5,6), `0000003000000` (6,5), `0000023000000`, `0000030200000`. That is 8 ends, so the class (mirror closed) cannot have an odd number. So 264 of 272 K = 5 ends, not 266 of 272, sit in rows with one long end. This does not touch the claim, only the "not carried by the lone 3" remark; it should be corrected (and it shows the lone 3 and its mirror are 2 of the 4).
Unchecked case: K >= 7 is empty in this class, so "K >= 6" is exactly "K = 6" (66 ends, 33 per side, one key). The "holds for K >= 6" wording is vacuous beyond K = 6, and the evidence for it is one K value in few distinct row shapes. The text says "(all K = 6)" but the title and claim should too.
The (4,7) class: 12 image keys at K >= 5, `removeVertex` None for one end. The author states this and does not claim it; fine. It also means "whole lone-3 key class (5,6)/(6,5)" in the title is one of two classes scanned, and the other does not obey the law at all. The title is narrow enough.

## New?
Grepped EXPERIMENTS.md for K0, E-146, E-173, lone 3. Predicted: E-146 (K0 = 5 from the single-3 table), E-135 (key of lone relation depends on {h,K}), E-141 (n = 13 refuted K0 = 4 at orbit level). Witness: E-173 (one row). Sizing: E-173 (2 674 440 rows), E-224/E-236 (shard plan). New: the full class scan (18 328 / 21 045 rows, K x key table) and the two-pass plan. The author's "Prior record" is accurate; "bears on: none" is right (no finding changes).

## Evidenced?
Mostly yes: table, examples, shard log and archive are committed and reproduce. Missing: the per-row shapes behind the K = 6 line (which of the 66 ends sit in rows with how many relations; 4 822 of the rows have 6 relations, so a typical K = 6 end belongs to a short core with no long end opposite); and the "refuted by" clause cites shard logs, not a count check that the 96 hits sum to 18 328 (I did that: ok).

## Scope
"key level only", n = 15, Coxeter key, are stated in title (key level), scope and claim. Title wording "its K >= 6 free ends delete to one n = 14 key" should say "K = 6 (no class member has K >= 7)".

## Required for acceptance
1. Correct "3 rows (6 ends)" and "266 of 272" to 4 rows (8 ends) and 264 of 272, naming the four rows.
2. State in title or claim that K >= 6 means K = 6 only, as the class has no K >= 7 end.
3. Add a note that shard files lack trailing newlines (concatenate with a separator or use `ends`).
