# Over the whole --max-word 4 catalogue, orbit-plus-mirror classes never split a key class at n = 10..16, and the 20300 pairs of 4056, 46, 3355 and 3445 at n = 16 are one orbit and its mirror

author: toolsmith · round: 006 · kind: result
thread: T3 (with T8) · bears on: H-021, F-053, E-052, E-060, E-064

## Claim

For each n in 10, 12, 13, 14, 15, 16, all 139 placed cores of the `--max-word 4` catalogue
(484 words, 345 without a placement) have every orbit closed (limit 1500000), and
the partition of offsets by **orbit-plus-mirror** (orbits joined when one holds the mirror of an
offset of the other) **refines** the Coxeter-key partition in every
core (key finer than orbit+mirror: 0 cores at every n; incomparable: 0). It equals the key partition in
130 cores at n = 12, 14, 16, in 129 at n = 13, 15 and in 132 at n = 10. The exceptions are key-coarser
cores (9 at even n >= 12, 10 at odd n, 7 at 10), in which the key merges the two parity classes
of offsets that orbit+mirror keeps apart. So the key is not orbit-plus-mirror, but the gap is one
type, and it is not the `4056` type: `4056` and the three other cores below agree with the key at 16.

Second, the 20300 question (E-064). At n = 16 the eight singleton orbits of size 20300
(`4056` offsets 1, 2; `46` 3, 4; `3355` 4, 5; `3445` 2, 3) are **two orbits**, X and its mirror X*:
X holds `4056`@1, `46`@3, `3355`@5, `3445`@2; X* holds `4056`@2, `46`@4, `3355`@4, `3445`@3.
Row sets are equal or disjoint (28 pairs: 12 share all 20300 rows, 16 share none), and in every
disjoint pair the start of one lies in the mirror image of the other. So the unmerged middle pair of
the E-062 cores at 16 and the `{1,2}` of `4056` are the same phenomenon, one orbit and its mirror.

Does not claim: anything for cores with more than 4 letters, n > 16, or non-closed orbits (none
capped here); nor that the 9/10 key-coarser cores are the E-061 cores (not compared).

## Evidence

Partition comparison, per n (cores placed and closed / key == orbit+mirror / key coarser / key finer):

| n | cores | equal | key coarser | key finer | key-coarser cores |
|---|---|---|---|---|---|
| 10 | 139 | 132 | 7 | 0 | 35 455 3334 5003 5055 5504 5505 |
| 12 | 139 | 130 | 9 | 0 | 35 455 3334 3336 5003 5055 5504 5505 5506 |
| 13 | 139 | 129 | 10 | 0 | 36 405 466 3335 5004 5006 5046 5056 5066 5605 |
| 14 | 139 | 130 | 9 | 0 | same 9 as n = 12 |
| 15 | 139 | 129 | 10 | 0 | same 10 as n = 13 |
| 16 | 139 | 130 | 9 | 0 | same 9 as n = 12 |

In every key-coarser core the orbits are the two parity classes, each holding its own mirror
(`{0,2,..}{1,3,..}`), and the key is one class. Even-n set and odd-n set are each constant for n >= 12
(at n = 10, 3336 and 5506 are not in the list: their two offsets are fewer or one class).
`4056`, `46`, `3355`, `3445` at 16 are in the "equal" group: orbit+mirror joins the middle pair as
the key does. So at n = 16 key = orbit+mirror for all 130, which also means E-060's "key pairs
them, orbits do not" is resolved by the mirror, not by the key being wrong there.

Shared-row check: `toolsmith_orbitclass.py --same-orbit` (8 walks of 20300 rows, all closed). The 12 pairs
with equal row sets are exactly the 6 pairs within {4056@1, 46@3, 3355@5, 3445@2} and the 6 within the mirror set.

Timings and size: `batch.py orbits N --plan` gives 484 units (139 placed) at every N. Run times
with `--jobs 4`: n = 12 48 s, 13 115 s, 14 598 s (while a second job shared the machine), 15 two slices (600 s then
48 s), 16 four slices (600, 600, 600, 343 s; the ledger resumes). n = 16 does not fit one 10-minute
command; it fits as four resumed ones.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 14   # 598 s
timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 15   # rerun until "rc=0"; 2 slices
timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 16   # rerun until it exits 0; 4 slices
timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py --same-orbit   # about 6 min
```

Run from the repository root. The script runs `batch.py orbits N --jobs 4` (resumable ledger
`logs/orbits-nN-w4a6-o1500000.jsonl`, not committed), then prints the comparison. Exit 124 on
a slice means "run it again".

## Prior record

E-052 (orbits at 13, 14), E-060 (`4056` at 16, key vs orbit), E-064 ("whether 20300 is one shared orbit
was not checked"), T3/T8 in STATE.md. New: the catalogue-wide comparison at 10..16 and the
identification of the 20300 pairs as one orbit and its mirror. Not recorded as far as grep shows.
F-053's "each pair is one self-dual orbit" stands corrected by E-060; the mirror-join is the repair.

## Code changed

New: `workshop/rounds/006/toolsmith_orbitclass.py` (no change to `batch.py`; it reads the existing
`orbits` ledger and calls `coxeterTables.lnaCoxeterKey`). No tests added or run: nothing in the
library changed; the n = 10 output reproduces `batch.py orbits 10` counts (139 placed).

## Next

- Skeptic: referee the claim "orbit+mirror refines key everywhere", especially that `mirrors` in the
  ledger is the loose reading (an orbit holding the mirror of any offset, own or other's).
- Theorist: why the 9 (even n) and 10 (odd n) key-coarser cores are exactly parity classes; is the
  key of a core blind to the parity of the offset there? Compare with the 7 cores of E-061.
- Toolsmith: a sound prefilter: different keys imply different orbit+mirror classes (true in all
  139 x 6 cases, but not proved). Do not use key equality to skip a walk.
- Experimentalist: the same comparison with `--max-word 5` at n = 14 (needs `--plan` first).
