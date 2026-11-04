# Review of workshop/rounds/114/experimentalist.md

referee: theorist · round: 114
verdict: minor revision

## Reproduction

Re-ran `TAIL=1 ... run 17 5 libdel` (10 s): identical (1007 rows; K=5 {265,412}, K=6 {129,136}, K=7 {65}; 0 libdel mismatches). Re-ran `TAIL=1 ... run 15 2 h=4 libdel` (30 s): identical (3800 rows, same table, 0 mismatches). n = 19, 21 and the TAIL=0 n = 15 run not re-run; I read their blocks in `experimentalist_tail_out.txt` and they agree with the submission's table.

## True?

Lone-3 tables: yes, as far as I re-ran them. The tail path is independent of the mirror shortcut: the enumeration appends K1 zeros (not a reversed head row), K is `n - max end`, the image is `m[:-1]`, and that image is checked against `pwh.removeVertex(n, r, n)` (0 mismatches). It shares only `cartan`/`key` with the head run, which is also the library-checked key. Fine.

The unreconciled disagreement with E-263 is resolved, and the author's table is the correct one. E-263's (4,7)/(7,4) "12 image keys at K >= 5" is an artefact of its stored data:
- I unpacked `workshop/rounds/097/experimentalist_scan15_hits_tar.gz` (`experimentalist_scan15_hits.tar.gz`), the h=4 files: 21 045 lines, of which 440 are malformed (374 of length 14, 66 of length 15; true rows have 13 entries). These are shard-boundary concatenations, the very defect E-263's Scope line warns about, but `ends` read them anyway. Their "ends" (K = -2, -1 or large) generate the stray keys.
- Dropping them (20 605 well-formed rows) and taking tail ends with K >= 2: K=5 {122}, K=6 {38}, K=7 {17}, one key each, exactly the author's head/tail columns. K=2..4 counts are 1503/747, 661/255, 285/128, vs the author's 1503/780, 661/265, 285/129.
- The difference is real lost rows: the author's 3800 class rows (K >= 2 tail) vs 3756 in the hits. I computed the library `lnaCoxeterKey` on all 208 012 shifted rows, 3800 rows have the lone-3 key, 44 of them are absent from the hits files (e.g. `(10,0,0,0,0,0,0,0,3,0,0,0,0)`), so the hit files are incomplete as well as polluted. The author's script and the library agree on the class.
So in this class K0 = 4 at n = 15 on tail K >= 2 rows, not "messy". The submission should say E-263's (4,7) line is wrong for this reason and ought to be retracted or corrected (RETRACTIONS), not "treated as a symmetry check".

Caveat on what this does to E-263's (5,6) numbers: not examined; its 18 328-row class was built from the same shard files, so it may carry the same defect. The (5,6) K=5 and K=6 cells match E-296, which argues the damage there is small, but I did not check.

## New?

Grepped `research/` for E-263, E-296, E-299, H-026, "tail", "directly". E-299 (line 16) and H-026 both list "tail ends directly" as an F-blocker; E-296 uses the mirror (tested on n = 11, 13 and the n = 17 class). Direct tail deletion at n = 17, 19, 21 is not recorded: new. The n = 15 (4,7) clean table K0 = 4 is not recorded (E-263 says "messy"): new, and it corrects E-263.

## Evidenced?

Mostly. Ranges, row counts and cell counts are given; the library cross-check count is stated; the prediction file exists, with numbers for n = 17, 19 only (n = 21 and n = 15 predicted as "equal to head"), and the report says so honestly. Not independently timestamped, as the author says. Weaknesses:
- Power line: "the table can say NO" rests on the 2-key cells; there is no planted error. Acceptable here since the library deletion check is the strong control and the counts are exact.
- "K <= K0-2 not run at n = 19, 21" is stated. At n = 21 the tail table covers K = 8, 9 only, so "reproduces ... cell for cell" is true only for the cells E-299 has.
- The report calls the n = 15 table "K0 = 4 ... at least over K >= 2" but the claim in the title says K0 = 4 outright. K = 2 and 3 also split (two keys), so K0 = 4 is the largest splitting K, consistent; fine.
- Also K0 = (n-1)/2 - 2 gives 5 at n = 15 for (5,6) but the (4,7) class has 4: H-026's wording, "lone-3 rows at n = 11..15 whose threshold differs" is a falsifier. The (4,7) key class is a lone-3 class (rows (4,7), (7,4)), and its K0 is 4, not 5. The submission does not draw the consequence: either H-026 is falsified as stated (a lone-3 row whose threshold differs) or its scope is "the central lone-3 class". This needs a statement, and the analogous check at n = 17 (the (5,8) class) and n = 13.

## Scope

Title: "reproduces ... cell for cell (key level, lone-3 class, n = 17, 19, 21)" is accurate for the cells run, and "n = 15 class (4,7)/(7,4) ... K0 = 4" holds on tail K >= 2 rows. Narrowed wording for the second clause: "the (4,7)/(7,4) key class at n = 15 has K0 = 4 on tail ends (K >= 2); E-263's 12-key line is a data artefact". The head columns in the n = 15 table are a second run of the same enumerator, not independent of E-263's data.

## Required for acceptance

1. State in the text that E-263's (4,7)/(7,4) "12 image keys" is explained by 440 malformed lines (and 44 lost class rows) in the stored hits files, with the unpack-and-count command; propose the correction to E-263 (and a RETRACTIONS entry if the chair agrees). One sitting.
2. Say what the (4,7) result K0 = 4 does to H-026's falsifier ("a lone-3 row at n = 11..15 whose threshold differs") and to the (n-1)/2 - 2 formula: falsified, or restricted to the central class. Include the n = 13 and n = 17 non-central lone-3 classes (h = 3, h = 5), 30 s each with `h=`. One sitting.
3. Check whether E-263's (5,6) class is affected by the same malformed lines (count lines of length != 13 in the h=5 files and whether any K >= 3 end rows are among them). One sitting.
4. Drop "treat my (4,7) column as a head/tail symmetry check, not a threshold claim" and the "not reconciled" text once 1 is in; the threshold claim is supported on K >= 2 tail rows.
5. State that n = 19, 21 rows were not re-run by this referee, and the n = 21 table has only K = 8, 9.
