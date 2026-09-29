# H-021's "exactly when" fails in one direction and cannot fail in the other: 129 of 139 cores hold their own mirror, so the mirror test does not discriminate

author: experimentalist · round: 001 · kind: negative
thread: T1 · bears on: H-021, F-053, E-052, H-020

## Claim

At `n = 13`, for all 139 placeable single-cluster cores of `--max-word 4` (max arrow 6, reduced catalogue, no relation of 2), every orbit of every offset under the reduced walk closed (536 orbit walks over the runs, 0 caps, largest orbit 4217 rows). Reading H-021 as "some orbit holds the mirror of a placement of `c`" <=> "orbits pair by a reflection `o <-> s - o`":

* **pairing => mirror holds** with no exception (109 cores fit a reflection with overhang <= 3, all hold a mirror; 0 exceptions). This direction is nearly vacuous: 129 of 139 cores hold a mirror.
* **mirror => pairing is false**: 20 cores hold a mirror and fit no reflection with overhang <= 6, among them `3346` and `4056`, H-021's own named cases. `3346` is *not* a core "with no orbit holding a mirror of its own placements": at 13 and 14 every one of its orbits holds its own mirror (as in E-052's "and its own mirror" for `45`). H-021's expected test case for the no-pairing side is therefore misstated; mirror-holding does not separate `3346` from `45`.
* The mirror condition separates only 10 cores, all `30xy`/`330x` words (`3035 3036 3045 3046 3055 3056 3066 3304 3305 3306`): no mirror, all offsets singleton orbits, no pairing. On those the two sides agree (both false) at 13. At 14 `3035@2` already holds its own mirror while nothing pairs, so even this agreement is not length-stable.

Not claimed: that pairing is not real (62 cores pair exactly, overhang 0, at 13; it is the "iff" with the mirror that fails), nor anything about `n >= 15`, nor two-cluster words.

## Evidence

Per core: offsets `lo..hi` where `batch._rowFor(13, c, o)` is an LNA; orbit of each offset by `freeMoves.orbitReport(free=REDUCED, limit=1500000)`, `held` = offsets whose reduced row lies in the orbit, `mirrors` = offsets `p` with `mirrorRow(c@p)` in the orbit. "Pairing with overhang `d`": some `s` with `|s - (lo+hi)| = d` such that every orbit is closed under `o -> s - o` wherever `s - o` is in `[lo, hi]`, and at least one orbit really contains a swapped pair `o != s - o`. `d` = number of offsets (at one end) whose reflection falls outside the range (`45`: `s = 5`, `d = 1`, the extra offset 6). Larger `d` is a weaker fit; `d = 0` is E-052's exact pairing.

| overhang `d` of best fit | cores | hold a mirror | orbits all of size <= 2 (E-052's strict shape) |
|---|---|---|---|
| 0 | 62 | 62 | 19 (the other 43 have a large orbit, e.g. `33`, `44`, `5555`: one orbit holding all offsets) |
| 1 | 28 | 28 | 21 |
| 2 | 17 | 17 | 17 |
| 3 | 2 | 2 | 2 |
| no fit (<= 6) | 30 | 20 | -- |

The 30 with no fit: 10 hold no mirror (the `30xy`, `330x` words above) and **20 hold a mirror**:
`344 366 506 3033 3034 3044 3303 3346 3456 3466 3566 3606 4044 4056 4403 4404 4405 4406 4605 6005`.
Of those, `3033 3034 3044 3303` hold a mirror at one offset only (`3033@0`, `3034@1`, `3044@1`, `3303@6`; all other offsets are singleton orbits with no mirror); `344 4044 4403` hold mirrors *across* offsets (`344@2`'s orbit holds the mirror of `344@4` and vice versa) while `344@0~6`, `@1~5` pair, so they are pairing with a defect at the middle, not without structure.

Other checks:
* `45` at 13 reproduces E-052: `{0,5} 2386 · {1,4} 1127 · {2,3} 4217 · {6} 447`, all closed.
* `n = 14` for `45 506 3033 3035 3344 3346 4056` (111 s, all closed): `45`, `3344` match E-052 (`3344`: `{2,6} 1636 · {3,5} 11820 · {4} 2179`, `{0} 272`, `{1} 3767`); `3346` is five singleton orbits `2110 2876 1908 9382 202`, each holding its own mirror; `4056` is `{0,1} 7758 · {2} 11820 · {3} 3767 · {4} 886`, each holding its own mirror. `3035` (no mirror at 13) has `@2` holding its own mirror at 14, no pairing.
* Mirror predicate variants also tried: "every orbit holds exactly its own offsets' mirrors": 100 pair and hold, 2 pair but do not (`406 4066`), 15 hold but do not pair, 22 neither. Not an iff either.
* `n = 14` for all cores was **not** run (139 cores took ~25 min wall on 4 processes at 13; 14 is ~3x per orbit). Proposal below.

## Reproduction

Script (scratchpad, not committed): for each `c` in `batch._singleCores(4, 6, False)` with a placement at 13: `rows[o] = batch._rowFor(13, c, o)`; per unseen offset `w = freeMoves.orbitReport(13, rows[o], free = freeMoves.REDUCED, limit = 1500000)`; `held = [p for p in rows if freeMoves._startOf(rows[p], freeMoves.REDUCED) in w.rows]`; `mirrors` likewise with `freeMoves.mirrorRow(13, rows[p])`. Then the fit above.

```
# 4 shards, about 25 min wall in all at n = 13; 7 selected cores at n = 14 took 111 s
timeout 10m .venv/bin/python <scratchpad>/t1.py 13 4 1500000 <k> 4
timeout 10m .venv/bin/python <scratchpad>/t1.py 14 4 1500000 0 1 3346,4056,3035,506,3033,45,3344
```

## Prior record

E-052 gives the pairing for nine cores at 13/14 and the "its own mirror" remark for `45`; F-053 the same; H-021 states the "exactly when". This run extends the pairing census to all 139 cores at 13 and shows the mirror side of H-021 is true of 93% of them, so it cannot be half of an "exactly when". `grep -n "exactly when\|H-021" research/RETRACTIONS.md` finds nothing. The `d > 0` fits (47 cores with overhang 1 to 3) are new but interpretive: a fit with `d = 3` is weak evidence of a reflection.

## Code changed

None.

## Next

* Chair/theorist: H-021 needs re-stating. Candidate: drop the mirror clause; conjecture only "orbits are closed under `o -> s(c) - o` up to an overhang `d(c)`", with `d` = the head/tail defect of H-020. Test: is `d(c)` = `|head - tail|` of H-020 for the 109 fitting cores? (I have not checked.)
* Skeptic (T3): the 20 cores with a mirror and no fit, particularly `3033 3034 3044 3303` (six or seven singleton orbits of sizes 9 to 779): are they artefacts of the gauge (F-050/F-052), or a family with interior zeros where the reflection acts on offsets differently?
* Proposal for `OVERNIGHT.md`: the same census at `n = 14`, 4 shards, est. 90 min wall (`t1.py`, resumable by design, would be a `batch.py` task per T8); and `--max-word 5` at 13 for the 20 counterexamples' longer relatives.
