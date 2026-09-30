# Experimentalist notebook (rewritten each round)

## What I now believe (after round 009)
- Mirror join done for 344/348/349 at n = 15..17 and 4046 at n = 14..16 (12 cells, all closed, `workshop/rounds/009/experimentalist_mirrorjoin_out.txt`): every equal-size singleton pair is one orbit plus its mirror, and orbit+mirror == key in all 12. Parity in these cores is orbit structure, plus the mirror.
- Key-coarser cores of E-064 (9 at n = 12: 35 455 3334 3336 5003 5055 5504 5505 5506; 10 at n = 13: 36 405 466 3335 5004 5006 5046 5056 5066 5605) are disjoint from the 7 of E-059 (344 366 4044 4403 4404 4405 4605). The 7 pair at even n and mirror-join at odd n; key-coarser cores are parity-class orbits each holding its own mirror. Not the same phenomenon.
- Earlier: k(34x) = x + 3, d = 0 (x = 4,5,7,8,9; n = 14..17); 346 one orbit; 45x no reflection; 4046 is a reflection (k = 11, d = 1); 5046/5056 translation at odd n; E-060's 4046@13 line did not reproduce.
- The size 20300 at n = 16 also occurs in 348 {2},{3} and 349 {1},{3}; shared rows with 4056's not tested.

## What I tried
- `workshop/rounds/009/experimentalist_mirrorjoin.py N WORD... [--limit L]` (wraps batch.orbitCensus and toolsmith classes). Timings: 348 n = 17 about 210 s, 349 n = 17 about 210 s, 4046 n = 16 about 170 s, 348 n = 16 about 260 s; running 5 jobs in parallel on 4 cores fit in 10 min.
- `toolsmith_orbitclass.py 12` and `13` (3-3.5 min each, ledgers in logs/, not committed).

## What I would do next
1. Key-coarser lists at n = 14, 15, 16 (is the even list stable across even n?). Needs `toolsmith_orbitclass.py 14` (maybe > 10 min: size with `--plan`; ledger resumes).
2. Is 20300 at 348/349 the same orbit X as in 4056/46/3355/3445 at 16 (`--same-orbit` style intersection)?
3. 34x at n = 18 and 44x remain unrun (E-068 open; chair said mirror-join first, now done).
4. Still open from 004: 3 failures at n = 17, 18; 139-core census n = 12 / 14 (n = 12 now done via orbitclass; fit not rerun).
- Watch: pair sums from classes with 3 members are not fits; sizes at n = 17 reach 120k, cap 300-400k.
