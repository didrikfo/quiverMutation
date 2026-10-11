# Equal-size singleton pairs of 344/348/349 (n = 15..17) and 4046 (n = 14..16) are each one orbit and its mirror; the key-coarser cores of E-066 at n = 12, 13 are disjoint from the 7 of E-061

author: experimentalist · round: 009 · kind: result
thread: T1/T3 · bears on: H-021', E-061, E-066, E-070, F-053

## Claim

In all 12 cells asked for (344, 348, 349 at n = 15, 16, 17; 4046 at n = 14, 15, 16) every orbit walk closed, and the orbit-plus-mirror partition equals the Coxeter-key partition. Each unmerged equal-size singleton pair is joined by the mirror: orbit X holds the mirror of the offset of orbit Y, so they are one orbit and its mirror, not two unrelated orbits of equal size. Second: the key-coarser cores of E-066 are 9 cores at n = 12 (`35 455 3334 3336 5003 5055 5504 5505 5506`) and 10 at n = 13 (`36 405 466 3335 5004 5006 5046 5056 5066 5605`); neither set meets the 7 cores of E-061 (`344 366 4044 4403 4404 4405 4605`), and the even and odd sets are disjoint from each other. Not claimed: that the mirror exhausts the odd-n defect in other cores (only 4 words were joined directly, plus the catalogue at n = 12, 13); nothing for n >= 18 or x >= 10.

## Evidence

Offsets and sizes (orbit size = rows in the reduced walk); `+mirror` is the join; all cells closed, limit 400000.

| n | core | orbits (size) | + mirror join | = key |
|---|---|---|---|---|
| 15 | 344 | {0,8}52 {1,7}43 {2}64 {3,5}51 {4}68 {6}64 | {2,6} joined, {4} centre | yes |
| 16 | 344 | 5 clean pairs | same | yes |
| 17 | 344 | {2}78 {8}78, {4}86 {6}86, {5}34 centre | {2,8}{4,6}{5} | yes |
| 15 | 348 | {0,4}8794 {1,3}4002 {2}4255 | same (centre) | yes |
| 16 | 348 | {0,5}77735 {1,4}19798 {2}20300 {3}20300 | {2,3} | yes |
| 17 | 348 | {2}9622 {4}9622 {3}9135 | {2,4}{3} | yes |
| 15 | 349 | {0,3}18416 {1,2}8693 | same | yes |
| 16 | 349 | {1}20300 {3}20300 {2}6860 | {1,3}{2} | yes |
| 17 | 349 | {1}20769 {4}20769 {2,3}20470 | {1,4} | yes |
| 14 | 4046 | {1}388 {2}388 {0,3}383 {4}11820 | {1,2} | yes |
| 15 | 4046 | {1}149 {3}149 {2}85 {0,4}81 {5}8794 | {1,3} | yes |
| 16 | 4046 | {1}530 {4}530 {2}534 {3}534 {0,5}521 {6}77735 | {1,4}{2,3} | yes |

Also 344, 4046 at n = 12 and 344, 348, 349, 4046 at n = 14 (pairs or mirror-joined, all = key). Every joined pair is also an equal-size pair, and at n = 16 the size 20300 of 348 {2},{3} and 349 {1},{3} is the same number as the 20300 of `4056`, `46`, `3355`, `3445` (E-066), a coincidence of size I did not test for shared rows.

Key-coarser vs the 7 (toolsmith's catalogue, `--max-word 4`, all 139 placed cores, all orbits closed): n = 12: orbit+mirror == key in 130, key coarser in 9 (list above), finer 0, incomparable 0. n = 13: 129, 10, 0, 0. In every key-coarser core the orbits are the parity classes `{0,2,..}{1,3,..}` (or a truncation) and each holds its own mirror, so the key merges them. For the 7: at n = 12 they pair cleanly (344 `{0,5}{1,4}{2,3}`), at n = 13 the singletons `{2}`, `{4}` hold each other's mirror (344 at 13: `{2}` mirrors `[4]`), and orbit+mirror == key. The 7 pair at n = 12 and mirror-join at n = 13, but the table shows mirror-joins at even n too (348@16, 349@16, 4046@14, 4046@16), so the "pair at even n, mirror-join at odd n" summary is NOT supported (chair, after referee); the key-coarser cores are parity-class orbits, never joined. Two different phenomena with a different list.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/009/experimentalist_mirrorjoin.py 17 344 348 349 --limit 400000   # about 7 min; n = 16 about 6 min; 4046 at 16 about 3 min
timeout 10m .venv/bin/python workshop/rounds/009/experimentalist_mirrorjoin.py 14 344 348 349 4046
timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 12 --jobs 4   # 3 min
timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 13 --jobs 4   # 3.5 min
```
Output of the first script, all cells: `workshop/rounds/009/experimentalist_mirrorjoin_out.txt` (n = 14 4046 line at the end).

## Prior record

E-066 has the mirror join for the catalogue at n = 10, 12..16 (and lists "key-coarser" counts 9/10, with the set not compared to E-061: its Limits). E-070 lists the singleton pair check for these cores as open; E-061 and E-070 say "equal size" only. New here: the direct join for these four cores at the named n, and the comparison of sets (n = 12, 13 names). The n = 12, 13 key-coarser lists were reproduced, not new (the counts 9/10 match E-066). The n = 14..16 key-coarser lists were not printed this round.

## Code changed

New `workshop/rounds/009/experimentalist_mirrorjoin.py` (uses `batch.orbitCensus` and `toolsmith_orbitclass.classes`). No library change, no tests run.

## Next

- theorist: characterise the key-coarser cores (parity-class orbits that hold their own mirror: all contain an orbit `{0,2,..}`); a pattern in the lists (`5xxx` with few letters, `35/36/455/466`) is visible but the lists differ by parity of n, so some dependence on n mod 2 (possibly in how the word's letters sum) is needed.
- experimentalist (next): print the key-coarser lists at n = 14, 15, 16 and check the even list repeats at 14 and 16 and the odd at 15; ask whether any key-coarser core has a pairing at another n.
- skeptic: is key-coarser a function of word and parity only? A referee reruns 348@16 (4 min).


## Chair note (round 009, after referee)
Title and claim 2 are scoped to n = 12, 13 (the n = 14..16 key-coarser lists were not printed). The direct joins at n = 14..16 duplicate E-066 rows; only n = 17 lies outside E-066's range. 344@16 sizes are in `experimentalist_mirrorjoin_out.txt`.
