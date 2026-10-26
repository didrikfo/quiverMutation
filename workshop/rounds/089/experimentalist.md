# |orbit(B_n)| = n 2^(n-5) holds at n = 13..16 as well (forward orbit, closed); A_n not in B_n's orbit (trivial: A_n has no move)

author: experimentalist · round: 089 · kind: result
thread: S-1 · bears on: E-218 (orbit-size fit), H-020 not touched
scope: n = 6..16 (E-218 covered 6..12; new: 13..16). Forward orbit of B_n = row 0..0 3 0 (kA_n/rad^(n-1) analogue with one 3-arrow end relation) under the same move set as E-218 (`fm.orbitReport(..., free=True)`, mirror-closed), all orbits closed (frontier emptied, limit 200000 not reached). Join test is membership of A_n or its mirror in the closed forward orbit of B_n; A_n has no move (single maximal relation rad^(n-1)), so non-joining is trivial (the content is the B_n orbit count). Says nothing about derived equivalence of A_n, B_n (only known at n = 9, 10).

## Response to referee

1. Mirror-closed column: added to the Evidence table (yes at every n = 6..16, from the script's `mirrorclosed` field).
2. Caveat moved into the title and first line of Scope/Claim: "A_n not in B_n's orbit (trivial: A_n has no move)"; the word "joined" is dropped.
3. n = 16 run: |orbit(B_16)| = 32768 = 16 * 2^11, closed, mirror-closed, A_16 not in it; 33.6 s. The prediction survived; claim extended to n = 16 and refutation now at n >= 17 (69632 predicted).
Reproduction notes (free=False, mirrorRow(A_13) = A_13) need no action.

## Claim

The E-218 fit holds three lengths further: |orbit(B_13)| = 3328 = 13 * 2^8 (prediction exactly met), and also 7168 at n = 14, 15360 at n = 15 and 32768 at n = 16 (the prediction, met exactly), matching n 2^(n-5). The orbit is closed and mirror-closed at every n = 6..16 and contains neither A_n nor its mirror, so A_n is not in the forward orbit of B_n or its mirror (trivially, A_n has no move). A_n has orbit 1 (no move), so non-joining is immediate: it is the B_n orbit count that carries content. Still a fit to 11 lengths, unproved. Refuted by: any n >= 17 with a different count (n = 17 predicted 69632, not run).

## Evidence

| n | orbit A_n | orbit B_n | predicted | match | A_n in orbit(B_n)+mirror | mirror-closed | time |
|---|---|---|---|---|---|---|---|
| 6 | 1 | 12 | 12 | yes | no | yes | 0.0 s |
| 7 | 1 | 28 | 28 | yes | no | yes | 0.0 s |
| 8 | 1 | 64 | 64 | yes | no | yes | 0.0 s |
| 9 | 1 | 144 | 144 | yes | no | yes | 0.0 s |
| 10 | 1 | 320 | 320 | yes | no | yes | 0.1 s |
| 11 | 1 | 704 | 704 | yes | no | yes | 0.4 s |
| 12 | 1 | 1536 | 1536 | yes | no | yes | 1.0 s |
| 13 | 1 | 3328 | 3328 | yes | no | yes | 2.7 s |
| 14 | 1 | 7168 | 7168 | yes | no | yes | 4.9 s |
| 15 | 1 | 15360 | 15360 | yes | no | yes | 14.5 s |
| 16 | 1 | 32768 | 32768 | yes | no | yes | 33.6 s |

Rows 6..12 reproduce E-218. Recurrence a_n = 2 a_{n-1} + 2^(n-5) also holds (3328 = 2*1536 + 256; 7168 = 2*3328 + 512; 15360 = 2*7168 + 1024). Output: `workshop/rounds/089/experimentalist_orb13_out.txt` (n <= 15), `..._orb16_out.txt` (n <= 16 rerun; same script, argument 16). Also 32768 = 2*15360 + 2048. No shard plan needed: the whole job takes about 25 s, nowhere near the 10 min limit.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/089/experimentalist_orb13.py 15   # ~25 s total; argument 16 gives ~70 s
```

## Prior record

E-218 (n = 6..12 fit, stated unproved, "forward orbit only" at n = 11, 12). grep of `research/` for n 2^(n-5) finds only E-218. Not new: the formula and the non-joining; new: n = 13, 14, 15 data points. A_n/B_n derived equivalence at n = 13 is not decided here (key equality only, as in E-218; not rechecked).

## Code changed

None in the package. New: `workshop/rounds/089/experimentalist_orb13.py` (+ output file). No tests.

## Next

Theorist: prove the recurrence (the orbit of B_n looks like 2^(n-5) copies indexed by n positions of a marker; a bijection would do). n = 17 (69632 predicted, ~2 min) can be run. Real open item remains whether A_n ~ B_n at n >= 11, which orbit counts cannot answer.
