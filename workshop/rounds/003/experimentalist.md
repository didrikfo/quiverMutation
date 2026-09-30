# The 7 "strict-mirror" cores pair at every even n from 12 to 18 and fail to pair at n = 13, 15, 17: the defect is a parity effect, not a property of n = 13

author: experimentalist · round: 003 · kind: result
thread: T1 · bears on: H-021, F-053, E-056

## Claim

At n = 14 all 7 survivors (`344 366 4044 4403 4404 4405 4605`) are reflection-pairing (fit d = 0 for 344 366 4403 4605, d = 1 for 4044 4404, d = 2 for 4405), hold no strict mirror, all orbits closed (26 orbit walks, 0 caps). So at n = 14 they are *not* mirror-without-reflection. Extending the same runs to n = 8..18: the 7 are strict-mirror-without-pairing at exactly n = 13, 15, 17 (7 of 7 each time) and pairing without strict mirror at n = 12, 14, 16, 18 (7 of 7 each time). Not claimed: that this parity holds for n > 18 or below 12 (n = 9..11 are mixed and have few offsets, so the fit is weak there); that the other 132 cores behave alike (unrun at n != 13, the full n = 14 census is parked); why parity matters (no mechanism offered); that the fit rule with d <= 6 is more than a descriptor.

## Evidence

Pairing status by n, the 7 cores, same rule as round 002 (`experimentalist_fit.py`; smallest d in 0..6). "Not placeable" cores at small n are dropped, hence fewer cores.

| n | cores | pair | strict mirror | no fit |
|---|---|---|---|---|
| 8 | 4 | 2 | 0 | 4044 4404 |
| 9 | 7 | 2 | 0 | 366 4044 4404 4405 4605 |
| 10 | 7 | 6 | 0 | 4405 |
| 11 | 7 | 6 | 0 | 4405 |
| 12 | 7 | 7 | 0 | none |
| 13 | 7 | 0 | 7 | all 7 |
| 14 | 7 | 7 | 0 | none |
| 15 | 7 | 0 | 7 | all 7 |
| 16 | 7 | 7 | 0 | none |
| 17 | 7 | 0 | 7 | all 7 |
| 18 | 7 | 7 | 0 | none |

Orbits of `344`: n = 14 `{0,7}47 {1,6}106 {2,5}114 {3,4}118` (four clean pairs, s = 7); n = 13 `{0,6}42 {1,5}34 {2}50 {3}19 {4}50` (the singletons 2 and 4 hold each other's mirror only); n = 15 `{0,8}52 {1,7}43 {2}64 {3,5}51 {4}68 {6}64` (again 2, 6 with equal sizes 64: singletons that would pair). Pattern: at odd n >= 13 the offsets that "should" be reflected split into singleton orbits of equal size (n = 13: 50 and 50; n = 15: 64 and 64), i.e. the reflection is present as an equality of orbit sizes but not as one merged orbit. At even n the same offsets merge into one orbit (`{2,5}114`). Orbit sizes for `344` at n = 12, 14, 16, 18: `{0,s}` 37, 47, 57 (n = 12, 14, 16); at odd n the `{0,s}` orbit is 42, 52 (n = 13, 15).

The n = 13 rows are re-derived from the committed round 002 census, not re-run.

## Reproduction

```
# each core at each n is 1 to 3 s; the whole sweep n = 8..18 is under 2 min
for n in 8 9 10 11 12 14 15 16 17 18; do for c in 344 366 4044 4403 4404 4405 4605; do timeout 10m .venv/bin/python workshop/rounds/002/experimentalist_census.py $n --cores $c --out logs/r3_${n}_$c.jsonl & done; wait; done
# (n = 13 rows come from workshop/rounds/002/experimentalist_census_n13.jsonl, filtered to the 7 cores)
.venv/bin/python workshop/rounds/003/experimentalist_table.py workshop/rounds/003/experimentalist_census_7cores_n8_18.jsonl   # seconds
.venv/bin/python workshop/rounds/002/experimentalist_fit.py workshop/rounds/003/experimentalist_census_7cores_n8_18.jsonl    # per-core lines, all n mixed
```
The census script has no `--plan` flag (the assignment mentioned one); sizing was by running one core: 3.4 s at n = 14. The saved file holds 74 rows (7 cores x n = 9..18, 4 at n = 8).

## Prior record

E-056 (the 7 at n = 13, "pairing with a defect" left open), T1 in STATE.md ("run the 7 at n = 14"). grep of `research/` for parity/odd n gives nothing relevant; nothing in RETRACTIONS. The answer to the T1 question is: not mirror-without-reflection at n = 14; and the n = 13 "defect" recurs at n = 15, 17, so E-056's 7-core list is a list of cores at odd n, not a family of cores that never pair.

## Code changed

Added `workshop/rounds/003/experimentalist_table.py` (per-n summary using round 002's fit script) and `experimentalist_census_7cores_n8_18.jsonl`. No library changes, no tests run or needed.

## Next

* Experimentalist (overnight, cheap): the full 139-core census at n = 14 will now tell whether the 30 no-fit and 13 mirror-without-fit cores also become pairing at even n, and n = 15 would show if the 62 exact cores lose pairing at odd n. If so, H-021 is a statement about parity of n and the "30 none" at 13 need rechecking at 12 and 14. Suggest adding n = 12 and 14 (not 15) to OVERNIGHT Menu 4 first: each is ~90 min.
* Theorist (conjecture, untested): the singletons 2 and 4 at n = 13 are the two half-orbits around a middle. A proof or a reason the middle offset splits at odd n >= 13 but not n = 9, 11 (where 4405 is the only failure) would restate H-021 with a parity clause.
* Skeptic: the fit with d <= 6 is generous; check the even-n pairs are exact (d = 0 for 4 of the 7) which is the strong evidence.
