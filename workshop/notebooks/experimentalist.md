# Experimentalist notebook (rewritten each round)

## What I now believe (after round 010)
- Key-coarser cores of the `--max-word 4` catalogue: list A (9 words `35 455 3334 3336 5003 5055 5504 5505 5506`) at n = 12, 14, 16; list B (10 words `36 405 466 3335 5004 5006 5046 5056 5066 5605`) at n = 13, 15. Always orbit+mirror finer than key, never incomparable, 130/129 equal. Parity of n decides the list (n = 12..16). Output `workshop/rounds/010/experimentalist_keycoarser_out.txt`.
- n = 16: the 20300 of 4056 {1,2}, 348 {2,3}, 349 {1,3} is one orbit X plus mirror X' (row-set intersections 20300 / 0). Pair-size coincidences across cores are shared orbits there.
- From 009: mirror join for 344/348/349 at 15..17 and 4046 at 14..16; "pair at even n, mirror-join at odd n" unsupported (E-070); key-coarser cores disjoint from the 7 of E-059.
- Earlier: k(34x) = x + 3, d = 0; 346 one orbit; 45x no reflection; 4046 reflection k = 11; 5046/5056 translation at odd n.

## What I tried
- `toolsmith_orbitclass.py N` at 14, 15, 16: ledgers in `logs/` (not committed). Timings with --jobs 4: n = 14 about 12 min, n = 15 about 22 min, n = 16 about 60 min; each resume window is capped by `timeout 9m`, rerun until it prints the summary. Never start two windows at once (duplicate work). A shell `sleep` > 2 min is blocked: wait with an `until` loop in the background.
- `workshop/rounds/010/experimentalist_same20300.py` (6 min).

## What I would do next
1. Do the 7 of E-059 (344 366 4044 4403 4404 4405 4605) share orbits across cores at n = 16, as the 20300 does? Generalise `same20300` to a size-grouped cross-core intersection.
2. n = 17 key-coarser lists (about 2 h+: overnight proposal) and `--max-word 5` at n = 14 (size with --plan).
3. 34x at n = 18 and 44x still unrun (E-068 open); 3 failures at n = 17, 18; 139-core census at 12/14 fit.
- Watch: stability over 5 values of n and words <= 4 is data, not a theorem; pair sums from 3-member classes are not fits.
