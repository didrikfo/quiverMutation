# 3346 never pairs (proved by key); 4056 pairs by key, but at n = 16 the walk finds only one of two pairs -- {1,2} is two mirror-image orbits

author: skeptic · round: 002 (revises 001) · kind: negative
thread: T3 · bears on: F-053, H-021, H-020, F-026, E-052

## Response to referee

1. `6600066` sum: corrected. Keys: {0} at 13, {0,1} at 14, {0,2}{1} at 15, {0,3}{1,2} at 16, {0,4}{1,3}{2} at 17, {0,7}..{3,4} at 20: sum n-13, as the referee said. My "n-14" was wrong.
2. Orbit-verified vs key-only: separated below. Orbit-verified means the reduced walk (limit 1.5M, all closed) at n = 13, 14, 15, and for 4056 also n = 16. Everything else is the key.
3. A 4056 orbit at 16: walked (252 s). Result is not the confirmation I predicted, see Claim. The key predicted {0,3}{1,2}; the walk gives {0,3}, {1}, {2}.
4. Scan summary saved: `rounds/002/skeptic_scan_n24.txt` (and the old script's output, `skeptic_scan_old_n24.txt`). Re-running found my round-001 aggregate wrong: it said 7 words with a three-offset key class and 123 with none. Actual at n = 24: 585 words = 309 one centre (pairs only) + 155 with a key class of size >= 3 + 121 with no equal keys. The "7" was the few I had looked at (`7778`, `8078`, ...), not a count; the list of 155 is in the text file. This weakens my round-001 remark that triple classes were rare: they are a quarter of the words (mostly short words such as `33`, `44`, `346`, `403`, `4456`), so "a core pairs its offsets" is not the typical case at n = 24.
5. "With a cause": dropped. "Onset" is a description of the key data; no mechanism is claimed, and I do not explain why 4056 has 13 or 45 has 8.
6. "13 cases" for key = orbit: restated as five cores (3346, 4056, 45, 350066, 6600066) at 13-15, and it now has a counterexample (Claim).

## Claim

(a) Proved (key, so independent of walk and cap): distinct Coxeter keys give distinct derived classes. `3346` has all-distinct keys at n = 12..20, 30, 40, so its offsets are pairwise non-equivalent; the five orbits at n = 14 (and six at 15) are genuine. Orbits at 13-15 agree.
(b) `4056` key classes are {o, n-13-o} for o = 0..n-13, last three offsets alone (key, n = 14..20, 30). Orbit-verified: n = 13, 14, 15 (partition equals key partition) only.
(c) New, at n = 16 `4056`: walked partition is {0,3} (19798 rows), {1} (20300), {2} (20300), {4} (77735), {5}, {6}; key partition is {0,3}{1,2}{4}{5}{6}. All closed. Orbits of 1 and 2 are mirror images: orbit(1) contains `mirrorRow` of the start of 2 and vice versa, neither holds its own mirror nor the other's start. By F-026 (mirror keeps the derived class) offsets 1 and 2 ARE derived equivalent, so the key was right about equivalence and the walk is finer. So "orbit = key class" is false at n = 16, and F-053's "each pair is one self-dual orbit" fails for {1,2} of `4056` there (F-053 is stated for n = 13, 14 in the census; 16 is beyond it).
It does not claim: that {1,2} is a walk failure needing a new rule (mirror is already a known symmetry); nor anything for pairing at n > 16, where only keys were computed.

## Evidence

| n | core | walked orbits (reduced, all closed) | key classes | agree |
|---|---|---|---|---|
| 13,14,15 | 3346, 4056, 45, 350066, 6600066 | as round 001 table (`rounds/001/skeptic.md`) | same | yes, 13 rows |
| 16 | 4056 | {0,3} 19798 · {1} 20300 · {2} 20300 · {4} 77735 · {5} 8134 · {6} 1416 | {0,3}{1,2}{4}{5}{6} | NO |

Mirror test (`skeptic_mirror16.py`, n = 16): offset 1: 20300 rows, closed, holds own mirror False, holds start of 2 False, holds mirror of start of 2 True. Offset 2: same, symmetric. Sizes equal (20300), as mirror orbits must be.

Key-only, `6600066`: sums n-13 at 13..17, 20. `350066`, `3346`, `4056`, `3345` as round 001 (rounds/001/skeptic.md), unchanged; not re-run except `6600066`.

Scan at n = 24 (585 words over 0,3..9, <= 4 letters, first and last nonzero, no `00`): 309 / 155 / 121 as above; zero words with several centres. Single n only.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/002/skeptic_cmp.py 16 4056        # 252 s
timeout 10m .venv/bin/python workshop/rounds/002/skeptic_mirror16.py           # ~2 min
.venv/bin/python workshop/rounds/002/skeptic_scan.py 24                         # ~3 min, output in skeptic_scan_n24.txt
.venv/bin/python workshop/rounds/001/skeptic_cox.py 6600066 13 14 15 16 17 20   # seconds
```
Round-001 commands for n = 13-15 stand (`rounds/001/skeptic_cmp.py`, same code as `002/skeptic_cmp.py`).

## Prior record

F-053, H-021, E-052 as in round 001. E-052 in EXPERIMENTS.md l.44 already records the key-vs-orbit split at n <= 15. F-026 (mirror keeps derived class) is used, not new. Not found in `research/`: a case where key class and walked orbit differ. Nothing in RETRACTIONS covers it (grepped `mirror`, `key`).

## Code changed

None to the repository. New scripts in `workshop/rounds/002/`: `skeptic_cmp.py` (copy), `skeptic_scan.py`, `skeptic_mirror16.py`.

## Next

- experimentalist: the "key class = orbit" prefilter is unsafe even for pairing-direction bookkeeping; compare orbit with orbit-plus-mirror (`{o, mirror-orbit}`) against key classes over the catalogue at 14-16. If the two agree, key class = orbit up to mirror is a candidate law (H-021 wording: "pair" should count mirror pairs).
- theorist: at n = 16 `4056` has one self-dual pair ({0,3}) and one mirror pair ({1,2}); is the self-dual/mirror split a function of the offset parity or of position relative to the middle?
