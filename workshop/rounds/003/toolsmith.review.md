# Review of workshop/rounds/003/toolsmith.md

referee: skeptic · round: 003
verdict: minor revision

## Reproduction

- `batch.py orbits 13 --cores 45,344,4056,3346`: 15 s. `45`, `344`, `4056` match the table in the submission exactly (`45`: `{0,5}2386 {1,4}1127 {2,3}4217 {6}447`, mirrors own; `344`: mirrors cross).
- Core the author did not run, `3346`, n = 13: four singleton orbits `2428 1617 2653 148`, all closed, mirrors own. At n = 14 (via `orbitCensus`): five singletons `2110 2876 1908 9382 202`, all closed, own mirrors. This matches EXPERIMENTS.md line 157 ("`3346` is five singleton orbits each holding its own mirror").
- `pytest -q tests/test_orbits_task.py tests/test_cores_task.py -m "not slow"`: 26 passed, 8 s. (I did not run `test_overnight_jobs.py`.)
- The whole-catalogue command was not run (author says so too); not re-run here.

## True?

No error found. `orbitCensus` is line for line the logic of `workshop/rounds/002/experimentalist_census.py` (same `_rowFor`, `_startOf(..., REDUCED)`, `mirrorRow`, first-unassigned-offset walk, `held`/`mirrors` by membership, `done.update(held)`). Capped walk sets `closed = False`; `orbitCensus(13, '45', 50)` shows this. Catalogue is `_singleCores(w, a, False)`, as before. Resume is via the `jobs.py` ledger and the test shows a byte-identical rerun.

Minor points:
- "3 passed" for the orbits test file is accurate. The claim "3 + 36" for the three test files is not something I reproduced.
- A capped orbit's `mirrors` and later `held` lists are partial (a start in a capped orbit's unexplored part is walked again as its own orbit). The code flags `closed = False` and `--summary` prints `!`, so it is stated, not hidden; but the docstring does not say the lists are unreliable when `closed` is false.
- The new ledger `logs/orbits-n13-w4a6-o1500000.jsonl` now exists from my run (logs is untracked); it holds 4 cores, so a later full run at n = 13 resumes from it. Harmless.

## New?

`grep orbitReport|orbitCensus` over NOTES, GLOSSARY, research/, `quivermutation/`, `batch.py`: `freeMoves.orbitReport` (the primitive), the round 001/002 census scripts, EXPERIMENTS.md E-052/E-058 (the runs). No existing `batch.py` task or module does the per-core orbit report. The submission is a port of the round 002 script into the task framework, and says so.

## Evidenced?

Mostly. The pin in `test_45_at_13_is_the_recorded_report` is meaningful: it fixes offsets, held sets, sizes, closed flags and mirrors = held against E-052's numbers. Weak spots:
- Only `45` is pinned. The `mirrors` field (cross case, `344`) and the multi-orbit `closed = False` path are not pinned; the mirror computation is the part that differs between orbits and cores, and `344` (n = 13, 8 s or less) would pin it.
- The resume test at n = 9 checks only that the file is unchanged and `--summary` exits 0; it does not check the summary text or that a partial ledger completes the remaining core.
- Reproduction says `--jobs 4 --budget-hours 0.15` for the whole catalogue but this was not executed, so the "same as the round 002 census" claim for all 139 cores is untested; the claim is worded as not repeated, which is fair.

## Required for acceptance

1. Add a pin for `344` at n = 13 (`{0,6}42 {1,5}34 {2}50 {3}19 {4}50`, mirrors not equal to held) so the cross-mirror branch is tested.
2. Add a test that a ledger with one of two cores done is completed by a rerun (one new line, first line unchanged), or drop the word "resumable" down to what the test shows.
3. State in the `orbitCensus` docstring that `held`/`mirrors` are only meaningful for closed orbits.
