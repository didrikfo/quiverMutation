# The key-coarser cores of the `--max-word 4` catalogue are the same 9 words at n = 12, 14, 16 and the same 10 at n = 13, 15; and the size-20300 orbits of `348`, `349` at n = 16 are the `4056` orbit and its mirror

author: experimentalist · round: 010 · kind: result
thread: T1/T3 · bears on: E-058, E-059, E-064, E-070, H-021'

## Claim

Over the 139 placed cores of the `--max-word 4` catalogue, every orbit closed (limit 1500000) at n = 14, 15, 16, and orbit-plus-mirror equals the Coxeter key in 130, 129, 130 cores; key coarser in 9, 10, 9; key finer or incomparable in none. The 9 key-coarser words at n = 14 and at n = 16 are exactly the 9 at n = 12 (`35 455 3334 3336 5003 5055 5504 5505 5506`); the 10 at n = 15 are exactly the 10 at n = 13 (`36 405 466 3335 5004 5006 5046 5056 5066 5605`). So the list is a function of the parity of n over n = 12..16 (n = 10 not rerun). In each, the orbits are the two parity classes `{0,2,..}{1,3,..}` (truncated for the longer words) and the key merges them; orbit+mirror leaves them unmerged. Second: at n = 16 the walks of `4056` {1}, `348` {2}, `349` {3} have identical row sets (20300 rows, shared 20300), and `4056` {2}, `348` {3}, `349` {1} likewise; the two sets are disjoint, and each contains the mirror of the other's start row. So the 20300 of `348`/`349` is the 4056 orbit X and its mirror X'. Not claimed: anything for n >= 17 or `--max-word 5`; whether the lists are stable for all even/odd n (the lists were stable over three even and two odd values only, and over words of length <= 4).

## Evidence

Key-coarser words (printed in full, with orbit/mirror/key partitions, in `workshop/rounds/010/experimentalist_keycoarser_out.txt`):

| n | placed+closed | orbit+mirror = key | key coarser | key finer / incomparable |
|---|---|---|---|---|
| 12 (E-070) | 139 | 130 | 9 (list A) | 0 / 0 |
| 13 (E-070) | 139 | 129 | 10 (list B) | 0 / 0 |
| 14 | 139 | 130 | 9 = list A | 0 / 0 |
| 15 | 139 | 129 | 10 = list B | 0 / 0 |
| 16 | 139 | 130 | 9 = list A | 0 / 0 |

List A = `35 455 3334 3336 5003 5055 5504 5505 5506`; list B = `36 405 466 3335 5004 5006 5046 5056 5066 5605`. Example partitions at n = 16: `35` orbits `{0,2,4,6,8}{1,3,5,7,9}`, key one block; `3336` `{0,2,4,6}{1,3,5}` | key `{0..6}`. The number of offsets grows by one per two units of n for each core, as the classes grow by period 2.

Same-orbit table (n = 16, size = rows of reduced walk, all closed): the six walks have size 20300; pair intersections are 20300 or 0:

| X (start rows) | X' = mirror of X |
|---|---|
| `4056`@1, `348`@2, `349`@3 | `4056`@2, `348`@3, `349`@1 |

Intersections between members of X are 20300 (a-start in b True); between X and X' 0, and the mirror of an X start lies in X'. Full output: `workshop/rounds/010/experimentalist_same20300_out.txt`.

Interpretation (a reading of the table, not a proof): at n = 16 the unmerged "middle pair" of `4056`, `46`, `3355`, `3445` (E-062/E-064), of `348` and of `349` (E-070) is one orbit X and its mirror, for all of them probably (`46 3355 3445` not rechecked here: E-064 identified their orbit with 4056's).

Cost: each of n = 14, 15, 16 took the ledger via `batch.py orbits` with `--jobs 4`: n = 14 about 12 min (two 9-min windows), n = 15 about 22 min, n = 16 about 60 min over 4 resumed windows of at most 9 min each (the ledger resumes; two overlapping windows at n = 15 duplicated work, no wrong rows because the reader keeps one record per unit: 484 distinct units at 14, 15, 16).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 14 --jobs 4   # repeat until it prints "n = 14: 139 cores"; ledger resumes
timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 15 --jobs 4   # about 3 windows
timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 16 --jobs 4   # about 7 windows (n = 16 units are slow)
timeout 10m .venv/bin/python workshop/rounds/010/experimentalist_same20300.py          # about 6 min
```
Do not start two windows of the same n at once (both recompute the same units).

## Prior record

E-064 gives the counts 9 (even n) and 10 (odd n) for n = 10, 12..16 but not the lists; E-070 printed the lists at 12 and 13 and left 14..16 open, and left the shared-rows test of 20300 open. E-064 had identified the 20300 pair of `4056` and the E-062 cores with one orbit and its mirror (E-070 recorded "same number, shared rows not tested" for `348`, `349`). New here: the lists at 14, 15, 16 and their equality with 12, 13, and the identification of `348`/`349` with the 4056 orbit by row intersection. No retraction concerns (RETRACTIONS not touched by this).

## Code changed

New `workshop/rounds/010/experimentalist_same20300.py`; no library change, no tests run. `toolsmith_orbitclass.py` run unchanged.

## Next

- theorist: why list A and list B (parity-class orbits that each hold their own mirror; two-letter and `5xxx` words); is the criterion "orbit parity classes at n, key = both" computable from the word?
- skeptic: is list A/B also stable at n = 17, 18 (orbitclass at 17 about 2+ h: overnight proposal) and at `--max-word 5` at n = 14 (needs `--plan`)?
- experimentalist: the 7 of E-059 against orbit X: do 344, 366, 4044... also share one orbit across cores at n = 16 (not asked).
