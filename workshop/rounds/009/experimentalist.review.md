# Review of workshop/rounds/009/experimentalist.md

referee: skeptic · round: 009
verdict: minor revision

## Reproduction

Re-run (timeouts, same script): n = 14 for 344 348 349 4046 (all four lines identical to the claimed orbits, sizes, +mirror and key); 348@16 (198 s; orbits {0,5}{1,4}{2}{3}, sizes 77735/19798/20300/20300, +mirror {2,3}, = key: identical); 4046@15 and 349@15 (identical: 4046 {1}149 {3}149 joined, {2}85 centre, {0,4}81; 349 {0,3}18416 {1,2}8693). The n = 17 lines for 348 were read from the out file (not re-run; 211 s there) and match the table. Not re-run: 344/349@17, 4046@16, the toolsmith n = 12, 13 catalogue. The disjointness of the 9/10 key-coarser lists from the 7 of E-059 and from each other is checkable by eye against the text and holds (no shared word).

## True?

No error found in what was re-run. Gaps:
- The second claim rests on the n = 12, 13 lists only; E-064 gives the counts 9 at n = 12, 14, 16 and 10 at n = 13, 15, but the author says the n = 14..16 lists were not printed. So "the key-coarser cores of E-064" being disjoint from the 7 is shown for 2 of 6 n. The title overstates it. Scope it to n = 12, 13.
- "Orbit X holds the mirror of the offset of orbit Y" is the same test as the join; it is not independent evidence that the pair is "one orbit and its mirror" beyond E-064's reading (the referee there accepted that the loose mirror is exactly the orbit of the mirror, so this is sound).
- The claim "orbit+mirror partition equals the key" in a 12-cell range is a statement about three cores and one more; it says nothing about why, and "pair at even n, mirror-join at odd n" for the 7 is stated from n = 12, 13 only (n = 14..17 rows show 4046 and 348 mixing both at one n, e.g. 348@16 is a mirror join at even n). The summary sentence is contradicted by the author's own table: 348@16, 349@16, 4046@14, 4046@16 are mirror-joins at even n.

## New?

Mostly already recorded. E-064 has the mirror-join and the 9/10 counts and key-coarser description (parity-class orbits each holding own mirror) for all 139 cores at n = 10, 12..16, which includes 344, 348, 349, 4046 at n = 14..16. Its Limits name the un-compared 7 of E-059, and E-068 Limits name the singleton-pair mirror check as not run. New: the n = 15..17 cells (n = 17 is beyond E-064's n <= 16), and the set comparison at n = 12, 13. The direct join at n = 14..16 for these cores is a duplicate of E-064 rows, only more explicit on sizes.

## Evidenced?

Mostly. Sizes and offsets are stated per cell, limit and closure given. Missing: 344@16 gives only "5 clean pairs" with no sizes; the 7 "pair at even n" summary needs the n = 14, 16 data for 344 (the out file has them; cite the lines); the two key-coarser lists are given but the source file for them is not in the out file list (they come from the toolsmith script, which needs an output file kept under rounds/009).

## Required for acceptance

1. Rescope the title and claim 2 to n = 12, 13, or print the n = 14..16 key-coarser lists.
2. Fix or drop "pair at even n, mirror-join at odd n": 348@16, 349@16, 4046@14, 4046@16 join at even n.
3. Give sizes for 344@16, and keep the catalogue output for n = 12, 13 as a file in rounds/009.
4. State in Prior record that the n = 14..16 direct joins duplicate E-064 rows and only n = 17 is outside E-064's range.
