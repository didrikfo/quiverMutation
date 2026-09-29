# 3346 and 4056 are not artefacts: the Coxeter key refutes their pairing outright, and 4056 is only F-053 below its onset length

author: skeptic · round: 001 · kind: negative
thread: T3 · bears on: F-053, H-021, H-020, F-050, F-052, E-052

## Claim

Neither exception to F-053 is an artefact of a cap, the walk or the gauge. (1) Every orbit
at issue closed (`stoppedBy == 'closed'`, largest 17905 rows against a 1.5M limit), so no cap
is involved. (2) The Coxeter key `coxeterTables.lnaCoxeterKey` -- a derived invariant, read off
the relation lengths, no walk and no sign gauge -- takes a *different value at every offset of
`3346`* at every `n` from 12 to 40. Two rows with different keys are not derived equivalent, so
the five orbits of `3346` at 14 are five genuinely different derived classes: not a failure of
the move set to join them. `3346` never pairs (checked to n = 40).
(3) `4056` is not "0 with 1, 2 alone" as a law; it is F-053's reflection with an onset: offsets
`0 .. n-13` pair with sum `n-13` and the last three offsets (`n-12 .. n-10`) are alone. At 14
that is {0,1}, 2, 3, 4; at 16 {0,3},{1,2},4,5,6; at 20 {0,7},{1,6},{2,5},{3,4}, 8, 9, 10.
Below n = 14 there is nothing to pair.
It does not claim: that the Coxeter partition equals the orbit partition in general (it does in
all 13 cases I compared, but the orbit is a subset of a key class, never the reverse, so for
*non-pairing* the key is a proof and for *pairing* it is only a necessary condition).

## Evidence

Orbit partition (reduced walk, all offsets, limit 1.5M) against the Coxeter-key partition:

| n | core | orbits (all closed) | key classes | agree |
|---|---|---|---|---|
| 13 | 3346 | {0}{1}{2}{3} | same | yes |
| 13 | 4056 | {0}{1}{2}{3} | same | yes |
| 13 | 45 | {0,5}{1,4}{2,3}{6} | same | yes |
| 14 | 3346 | five singletons (2110, 2876, 1908, 9382, 202 rows) | same | yes |
| 14 | 4056 | {0,1} 7758 · {2} 11820 · {3} 3767 · {4} 886 | same | yes |
| 14 | 45 | {0,6}{1,5}{2,4}{3}{7} | same | yes |
| 15 | 3346 | six singletons (up to 17905 rows) | same | yes |
| 15 | 4056 | {0,2} 4002 · {1}{3}{4}{5} | same | yes |
| 13 | 350066, 6600066 | {0}{1}; {0} | same | yes |

Key partitions only (no walk, seconds):

| core | pairing by key |
|---|---|
| 3346 | none, n = 12..20, 30, 40 |
| 4056 | pairs sum n-13, tail of 3 alone, n = 14..20, 30 |
| 350066 (an H-020 failure) | nothing to 16; {4,5} at 17; {4,6} at 18; {4,7}{5,6} at 19; {4,8}{5,7} at 20 |
| 6600066 (an H-020 failure) | {0} at 13, {0,1} at 14, then pairs sum n-14 |
| 3345 (control) | pairs sum n-9 at 12, 13, 14, 16 |

So the two H-020 failures I could name (`350066`, `6600066`) are the same phenomenon: the
slide at 13 is one or two offsets because the reflection has not yet got room. It is the
"three offsets or more" reading of H-020's second amendment, now with a cause: a core has an
onset length below which no two offsets share a class. I did not reach the other four; the
census ledgers are not in the repository, only their two examples in E-051.

Scan of every word over digits `0,3..9`, at most 4 letters, first and last nonzero, no `00`
(585 words at n = 24): 309 have exactly one pairing centre, 7 (e.g. `7778`, `8078`) have a key
class of three offsets, and 123 have no equal keys at all at n = 24 (`3346` among them,
also `304..309`, `3033..3099`, `3456`, `3566`, ...). Long list: `rounds/001/skeptic_scan.py`.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/001/skeptic_cmp.py 14 3346 4056 45   # 34 s + 16 s + 44 s
timeout 10m .venv/bin/python workshop/rounds/001/skeptic_cmp.py 15 3346 4056       # ~100 s
.venv/bin/python workshop/rounds/001/skeptic_cox.py 3346 12 14 20 40               # seconds
timeout 9m .venv/bin/python workshop/rounds/001/skeptic_scan.py 24                 # ~3 min
```

## Prior record

E-052 lists both as exceptions and F-053/H-021 leave them open. Nothing in `research/` uses
the Coxeter key to partition the offsets of a core (grep of `Coxeter key` in research/: used as a search gate and as a verifyMove check (E-033), never to
partition offsets -- I found no such use, but the grep is by phrase); F-050 (sign gauge) concerns the mutation
procedure, and the rule-table walk does not use it. F-052 (one-way free move) is why the
reduced walk was used here. Not in RETRACTIONS.

## Code changed

None to the repository. Three throwaway scripts saved as `rounds/001/skeptic_*.py`.

## Next

- experimentalist (T1): the Coxeter key is a cheap prefilter and a proof of separation.
  Before walking a core, partition offsets by key; walk only within key classes. Orbit =
  key class held in all 13 comparisons; test it across the whole `--max-word 4` catalogue at 13
  and 14 (the interesting case is a key class the walk splits, which would be a missing rule
  or a genuine non-equivalence with the same Coxeter polynomial).
- theorist (T2): the centre is `n - k(c)` with `k` a property of the core (`4056`: 13, `45`: 8,
  `3344`: 6, `556`: 10). Data for a formula: `skeptic_scan.py`. The 123 cores with no pair, and
  the 7 with triples, are the ones that would break a formula.
