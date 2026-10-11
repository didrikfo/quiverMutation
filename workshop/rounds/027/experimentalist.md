# Conjecture W (E-109) has 0 mismatches on 32 000 fresh out-degree 2 rows (n = 8 classes 0, 2, 3, n = 9 class 0 prefix, parallel-arrow rows included), but classes 2, 3 and n = 9 contain no out-degree 2 reject, so only class 0 tests it

author: experimentalist · round: 027 · kind: result
thread: T5 · bears on: E-109, E-110, E-111

## Claim

On four capped guarded BFS walks (540 s each, run 4 in parallel, so counts are load-dependent lower bounds on coverage) every gate-admitted
out-degree 2 row has W == (J != 0): 0 mismatches over 32 132 rows, 388 + 116 + 209 + 422 = 1 135 of them with the two out-arrows of v
parallel (same head), which E-109 skipped. Not claimed: that this tests the converse (J != 0 => W) outside n = 8 class 0: the only
out-degree 2 rejects are the 61 in class 0 (W true on all 61). Classes 2, 3 and n = 9 class 0 have 0 out-degree 2 rejects, so there
W == False == K is an agreement of negatives (no positive control). No parallel-arrow row, at any out-degree, was a reject (0 of 1 435
rows with parallel out-arrows, and the rows with parallel arrows elsewhere in the quiver are inside the other counts). Also new: out-degree >= 3
(3 494 rows, 300 parallel) has 0 rejects in every walk, and the "1 reject outside the long square" rows are only out-degree 1.

## Evidence

Key = (out-degree, parallel out-arrows at v, W, J != 0); rows are (algebra, v) pairs, gate-admitted, repeated across algebras.

| walk (540 s) | algebras expanded / seen | out 2 rows | of them parallel-out | W & J!=0 | W only | J!=0 only | out 1 J!=0 | out >= 3 rows | out >= 3 J!=0 |
|---|---|---|---|---|---|---|---|---|---|
| n=8 c0 (2 algebras in class) | 5 145 / 13 086 | 6 211 | 388 | 61 | 0 | 0 | 4 of 14 609 | 1 031 | 0 |
| n=8 c2 (18) | 9 093 / 22 638 | 9 989 | 116 | 0 | 0 | 0 | 0 of 25 993 | 1 071 | 0 |
| n=8 c3 (20) | 9 705 / 23 824 | 10 818 | 209 | 0 | 0 | 0 | 6 of 28 429 | 974 | 0 |
| n=9 c0 (2) | 4 553 / 10 555 | 5 114 | 422 | 0 | 0 | 0 | 110 of 15 945 | 418 | 0 |

(Only the first number of "expanded" is the coverage; "seen" is the distinct canonical keys found, the frontier left unexpanded.) n = 8 c0
shows 61 rejects in 5 145 algebras, against 42 distinct and 55 rows in E-108/round 026 on fewer algebras: rate about 1.2 % of out-2 rows,
again prefix dependent. The 4 + 6 + 110 out-degree 1 J != 0 rows are E-105's long-square family (not re-classified here).

Class sizes (from `toolsmith_rejwalk.py --plan`): n = 8 has 11 Coxeter classes, class 0 has 2 algebras, 1: 8, 2: 18, 3: 20; n = 9 has 19, class 0 has 2.
Coverage: n = 8 c2/c3 reached BFS depth with 22-24 k canonical keys found; n = 9 c0 10.5 k (E-107's 500 s walk saw 8 231). A 540 s cap
is a prefix; none of the walks closed. n = 9 c0 is the toolsmith's multi-night job (ratio 2.5 per level): not run here.

Weakness of the parallel test: W was evaluated through `relationsFrom` (uses `arrowRels`, keyed arrows), but since no parallel row rejects and
no relation of the kind W looks for appears with both out-arrows parallel (W False on all 1 135), this exercises the code path, not a positive.
A hand-built positive control (relation `p1 b1 = p2 b1`, parallel b1, b2 killing both) is the right test and is not written.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/027/experimentalist_w.py 8 2 540    # 9 min; args: n class budget_sec [max_exp]
timeout 10m .venv/bin/python workshop/rounds/027/experimentalist_w.py 8 3 540
timeout 10m .venv/bin/python workshop/rounds/027/experimentalist_w.py 9 0 540
timeout 10m .venv/bin/python workshop/rounds/027/experimentalist_w.py 8 0 540
```
Counts change with load; pass a 5th argument (max expansions, deterministic) for an exact rerun. Timing: about 40 algebras/s at n = 8 solo.

## Prior record

E-109 (W, 0 mismatches on 17 802 rows; limits: no n = 9, classes 2-3, out >= 3, parallel rows skipped), E-110, E-107 (n = 9 c0 prefix no
D' rejects), E-111 (rejwalk; n = 9 sized only). This submission fills the listed gaps in E-109's limits with null results: no refutation, and
no new positive. Nothing in RETRACTIONS touched.

## Code changed

New `workshop/rounds/027/experimentalist_w.py` (walk + inline W + kerdim; reuses `rounds/023/scholar_longsquare.py`). No library changes, no tests.

## Next

- toolsmith/overnight proposal: n = 8 c1 out-2 W check solo with `--max-exp`, then n = 9 c0 via `toolsmith_rejwalk.py ... --budget-hours 7 --ckpt`
  with W added to its tally (about 10 lines), to see whether an out-degree 2 reject ever appears at n = 9 (E-107 says none in 8 231 algebras; ours 10.5 k).
- theorist: the "cancels" branch of W is still unexercised; a hand-built algebra with W true only through a sum would test it.
- experimentalist: hand-built parallel positive control; classes 4..10 at n = 8 (sizes 26-266) for any out-degree 2 reject.
