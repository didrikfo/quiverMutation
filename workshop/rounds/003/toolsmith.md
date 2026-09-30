# `batch.py orbits` reproduces E-052's orbit report for `45` at n = 13 as a resumable, tested task

author: toolsmith · round: 003 · kind: tool
thread: T8 · bears on: H-021, F-053, E-052, E-058

## Claim

`.venv/bin/python batch.py orbits N [--cores ...]` now reports, per core, the
reduced-walk orbits of its placements at length N (held offsets, size, closed
or capped, and the offsets whose mirror each orbit holds), one core per ledger
line under `logs/`, resumable, with `--summary` printing a pairing/mirror
table. For `45` at n = 13 it gives `{0,5} 2386`, `{1,4} 1127`, `{2,3} 4217`,
`{6} 447`, each holding exactly its own offsets' mirrors, as recorded in
round 002; a test pins this (8 s, not slow). Not claimed: any new mathematics;
that the Coxeter key can serve as a prefilter (it cannot alone, E-058; none was
added); that the pairing summary decides H-021's overhang fit (it lists only
the sums `a+b` of two-offset orbits; the fit stays in `experimentalist_fit.py`).

## Evidence

Reuse audit of `workshop/rounds/002/experimentalist_census.py`: its logic is
`_rowFor` for each offset, reduced starts via `freeMoves._startOf`, mirrors via
`freeMoves.mirrorRow` (same as `batch._mirror`), an `orbitReport` walk from the
first unassigned offset, and `held`/`mirrors` by membership. I found no
soundness problem: the walk uses REDUCED, a start's own offset is always held
so the loop terminates, offsets are removed once held, capped walks are flagged
`closed = False`, and the catalogue is `_singleCores(w, a, False)` (no `2`
words, as the reduced walk requires). I moved that body into `batch.orbitCensus`
(about 30 lines), because `_rowFor` lives in `batch.py`; the census script's
other parts (own shard/JSONL resume) are replaced by the ledger of `jobs.py`.
The old script is untouched and its output format differs (ledger wraps each
record as `{unit, at, result}`; the fields inside are the same).

Checks run:

| check | result |
|---|---|
| `45`, n = 13 | `{0,5} 2386 {1,4} 1127 {2,3} 4217 {6} 447`, all closed, mirrors = held |
| `344`, n = 13 | `{0,6}42 {1,5}34 {2}50 {3}19 {4}50`; mirrors "cross" (orbit `{2}` holds the mirror of 4; round 002's description) |
| `4056`, n = 13 | four singletons `493 2116 2386 224`, no pairs, own mirrors |
| ledger resume | second run leaves the file byte-identical (test, n = 9) |

Timing: `45` alone 8 s; `45,344,4056` together 12 s (one process). The
round 002 census of all 139 cores at n = 13 took about 4 procs x minutes; the
task takes `--jobs` for the same.

## Reproduction

```
timeout 10m .venv/bin/python batch.py orbits 13 --cores 45,344,4056     # 12 s
.venv/bin/python batch.py orbits 13 --cores 45,344,4056 --summary       # instant
.venv/bin/python batch.py orbits 13 --jobs 4 --budget-hours 0.15        # whole catalogue, resumable (exit 2 = out of budget; rerun)
timeout 10m .venv/bin/python -m pytest -q tests/test_orbits_task.py -m "not slow"   # 3 passed, 8 s
```

The whole-catalogue command was not run here (round 002 recorded 139 cores at
n = 13 in the census; not repeated). For n = 14 (T9, overnight) the ledger name
carries n, so shards are `--cores` lists into the same ledger only if run one
at a time (the ledger is not safe for concurrent processes; use `--jobs`).
The default `--orbit-limit 1500000` matches the round 002 census; the ledger
name includes it.

## Prior record

E-052 (orbit pairing of `45`, first census), E-056 (strict/loose readings
census at n = 13, the source of the script), E-058 (key class is not orbit:
`4056` at n = 16). This adds no finding; it makes the census a task.
`research/` has no other orbit-report tool (grepped `orbitReport`, `census`
in `batch.py`/`README`: only the round 002 script).

## Code changed

- `batch.py`: `orbitCensus(length, word, limit)`, `OrbitsTask` (registered in
  `TASKS`), one line in the module docstring. No existing behaviour touched.
- `tests/test_orbits_task.py` (new): pins `45` at 13; no-placement returns
  `None`; task resume and `--summary` at n = 9.
- Ran `tests/test_orbits_task.py`, `tests/test_cores_task.py`,
  `tests/test_overnight_jobs.py` (`-m "not slow"`): all pass (3 + 36).
- Not changed: `research/`, `STATE.md`, the round 002 scripts.

## Next

- experimentalist: use `batch.py orbits 14 --jobs 4` for the 7 survivors
  (`--cores 344,366,4044,4403,4404,4405,4605`) instead of the old script.
- Prefilter: a key prefilter would only be safe as "skip the walk when
  orbit-plus-mirror is already known equal"; that needs a key-versus-orbit
  table over the catalogue at 14..16 (skeptic's request to experimentalist);
  I did not build it because no experiment needs it yet.
- toolsmith (later): `summarise` could print the reflection fit if T1's
  restatement fixes a definition.
