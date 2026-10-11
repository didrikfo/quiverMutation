# Out-degree-2 rejects with no long square recur at n = 8 class 1 but not (in these caps) at n = 8 class 3 or n = 9 class 0

author: experimentalist · round: 025 · kind: result
thread: T5 · bears on: E-102, E-105, H-015

## Claim

Using `workshop/rounds/023/scholar_longsquare.py` (guided walk, gate-admitted steps, wall-clock cap 500 s, four runs in parallel on 4 cores),
the J != 0 rows with out-degree 2 and no long square are: 48 at n = 8 class 0 (re-run; referee had 42 in 200 s), 15 at n = 8 class 1, 0 at n = 8 class 3,
0 at n = 9 class 0. All counts are cap-dependent lower bounds on a prefix of the walk, so the zeros are NOT evidence of absence: n = 9 covers only 8 231 algebras and
n = 8 c3 only 18 269 algebras with just 6 J != 0 rows in total (none out-degree 2).
Also new: n = 8 class 1 has 2 rejects with out-degree 1 and NO long square (a second way to break the iff), and at n = 8 c0 and c1 the out-degree-2 rejects are the majority of J != 0 rows.
Every J != 0 row in every run has tiltingPlus False; every J = 0 row has tiltingPlus True; `mono` is False in all rows. The J != 0 <=> tiltingPlus False link is intact; only the "long square" shape fails.
Not claimed: these are D-type (not classified); that they are reachable only at n >= 8 (n <= 7 data are from E-102/023 only); anything about completeness.

## Evidence

Coverage per run (algebras = distinct algebras visited by the capped walk; rows = (parent, v) gate-admitted steps, tallied by (J, outdeg, longsq)):

| run | algebras | rows total | J != 0 total | J != 0, out 1, longsq | J != 0, out 1, no longsq | J != 0, out 2, no longsq | J = 0 out 1 / out 2 |
|---|---|---|---|---|---|---|---|
| n=8 c0, 500 s | 10 478 | 17 565 | 50 | 2 | 0 | 48 | 11 882 / 5 633 |
| n=8 c1, 500 s | 16 976 | 27 839 | 21 | 4 | 2 | 15 | 19 042 / 8 776 |
| n=8 c3, 500 s | 18 269 | 31 192 | 6 | 6 | 0 | 0 | 22 573 / 8 613 |
| n=9 c0, 500 s | 8 231 | 16 557 | 85 | 85 | 0 | 0 | 12 277 / 4 195 |

(`rows total` summed from the tallies.) Class 2 at n = 8 and class 0 beyond the cap not run. Rates of J != 0 per row: 0.28 % (c0), 0.08 % (c1), 0.02 % (c3), 0.51 % (n=9 c0).
Reading: at n = 9 class 0 the 85 J != 0 rows all have out-degree 1 and a long square, as at n <= 7; at n = 8 c3 likewise (6 rows). The out-degree-2 rejects appear in the two
n = 8 classes where J != 0 is dominated by them. With 8 231 algebras at n = 9 (vs 10 478 at n = 8 c0, where the first out-degree-2 rejects appeared within the first 9 000), the n = 9 zero is weak:
the rows may simply come from a different part of the walk. No claim that n = 9 differs.
Caveat on timing: the four runs shared cores, so each covers fewer algebras than a solo run would; counts differ from the referee's 200 s solo run.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/023/scholar_longsquare.py 8 --class 0 --budget-sec 500   # 500 s (+ startup), 4 runs parallel
timeout 10m .venv/bin/python workshop/rounds/023/scholar_longsquare.py 8 --class 1 --budget-sec 500
timeout 10m .venv/bin/python workshop/rounds/023/scholar_longsquare.py 8 --class 3 --budget-sec 500
timeout 10m .venv/bin/python workshop/rounds/023/scholar_longsquare.py 9 --class 0 --budget-sec 500
```
Outputs: `workshop/rounds/025/experimentalist_n8_c0.txt`, `_n8_c1.txt`, `_n8_c3.txt`, `_n9_c0.txt`.

## Prior record

E-102 (n = 5..7 iff), E-105 and STATE T5 (42 out-degree-2 rejects at n = 8 c0, hasLongSquare is a one-out-arrow test). Not in RETRACTIONS. New: recurrence in n = 8 c1 (15), the
2 out-degree-1 no-long-square rejects at n = 8 c1 (a case that `hasLongSquare` misses even with one out arrow; possibly a presentation artefact like scholar's G; not inspected), and the null at n = 8 c3 / n = 9 c0 within caps.

## Code changed

None.

## Next

- Toolsmith: make the script print one example parent (arrows, rels, v) per tally key, so the 2 out-1 no-longsq rejects at n = 8 c1 can be classified (G-type presentation or genuinely new).
- Overnight proposal: n = 9 classes 0-1 with `--budget-sec 3000` solo (not parallel) and n = 8 c2, to see whether out-degree-2 rejects appear at n = 9 at all.
