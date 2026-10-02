# The 61 n = 8 class-0 D' rejects are admitted exactly because the gate tests single paths while the kernel element is a two-term sum; the D' shape is selective (61 of 84 rows with the loose shape) and occurs only in class 0

author: skeptic · round: 029 · kind: result
thread: T5 · bears on: E-105, E-110, E-111

## Claim

(1) Admission. The gate (`procedure.isMutable`) refuses v iff some nonzero single path p into v has p*beta = 0 for every out-arrow beta.
A reject is, by definition, J != 0 with no single path in J, so "admitted for the stated reason" is true, and it is close to a tautology:
the gate is the "path" version of J != 0, the reject is the "combination" version. The non-tautological content is the witness pattern
below (rejecting rows have paths witnessed by one arrow only, either one) and the shape. Not claimed: a proof of why only n = 8.
(2) Analogues. Out-degree 1 admitted rejects (the long-square family, E-097, E-103) are the same mechanism at out-degree 1 (x = p - q killed by the single arrow) and
occur at n = 6, 7, 8. At out-degree 2 the rejects occur only in n = 8 class 0 (61), in no other class or n in my sample.
(3) Null. The loose D' shape L (a two-term relation commuting into one out-arrow, x not in I, plus a monomial relation ending in the other out-arrow at v) is
NOT sufficient: it also appears at n = 6 c0 (5 rows), n = 7 c0 (20) and n = 8 c0 (23), all accepted. Adding the W condition (both terms die into the other arrow) separates perfectly.

## Evidence

Same guarded walk as E-111 (`scholar_longsquare.py` prefix, key-preserving, cycles skipped), capped by number of expanded algebras (counts are prefix-dependent).
Row = gate-admitted (algebra, v). L = loose shape; W = E-107 rule; K = J != 0 (kerdim). 53 162 algebras expanded in total.

| n, class (size) | algebras expanded | out-2 rows | L true | L true, K false | W & K | out-1 K |
|---|---|---|---|---|---|---|
| 8, c0 (2) | 5 145 | 6 211 | 84 | 23 | 61 | 4 |
| 8, c1..c3 (8, 18, 20) | 7 500 | see files | 0 | 0 | 0 | 0 |
| 8, c4..c10 | 12 500 | see files | 0 | 0 | 0 | 0 |
| 7, c0 (8) | 2 522 | 2 229 | 20 | 20 | 0 | 26 |
| 7, c1..c5 | 15 000 | see files | 0 | 0 | 0 | 0 |
| 6, c0 (2) | 3 000 | 3 008 | 5 | 5 | 0 | 15 |
| 6, c1..c3 | 9 000 | see files | 0 | 0 | 0 | 0 |

(Rows marked "see files": every key line has K false, L false; exact counts in the .txt files; out-2 row totals were not summed.)
Selectivity at n = 8 c0: of 6 211 out-2 rows, 84 have L (1.4 %); of those 61 reject (73 %), 23 accept; W & K 61, W & not K 0, not W & K 0.
So the D' shape is selective (1 % of rows, 0 rejects outside it), but L alone is not a criterion: the 23 L-true accepts (and 25 at n = 6, 7)
are the control the shape needs, and W's extra clause is what discriminates. Same-n comparison: n = 6, 7 c0 have L with 0 rejects.
Witness pattern at the 61 rejecting rows, over all nonzero paths into v (arrow 0 / arrow 1 / both with p*arrow nonzero): 172 / 85 / 63; no path has an empty set (that is the gate's admission,
re-derived independently of `mutationIsPossibleAtVertex` only in the sense that the same `isInIdeal` is used: not independent code).
Reproduces E-111 exactly at c0 (61 rejects, 4 out-1 rejects, 1 031 out >= 3 rows with none).

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/029/skeptic_dprime.py 8 0 0 5145      # 462 s under load (6 jobs in parallel); ~170 s solo
timeout 10m .venv/bin/python -u workshop/rounds/029/skeptic_dprime.py 8 1 3 2500       # also 4 6 2000, 7 10 1500; n = 7 and 6: `7 0 30 3000 100`
```
Output files `workshop/rounds/029/skeptic_dprime_n*.txt`. Counts other than the 61 depend on the cap.

## Prior record

E-105 (61 D' rejects; "why n = 8 and not n <= 7" open), E-106/E-108 (shape necessary not sufficient: 767 / 42), E-107, E-110 (J defined by the parent alone, so "reject => W" not derivable
from step 7), E-111 (0 out-2 rejects at c2, c3, n = 9 c0 prefix). New here: the loose D' shape exists, unrejected, at n = 6 and 7 class 0 (so the n <= 7 absence is not an absence of the shape),
classes 1, 4..10 at n = 8 and all classes at n = 6, 7 have no out-2 reject (rediscovery of E-111 for c1, c2, c3 and n <= 7, new for c4..c10 within the caps),
and the class-0-only count. Nothing in RETRACTIONS touched.

## Code changed

None in the library. New `workshop/rounds/029/skeptic_dprime.py`; no tests.

## Next

- Weak spots: the L column is my definition (loose D'), not E-105's; caps 1 500..5 145 algebras so classes 4..10 are prefixes, c1 and c2 at n = 8 not walked to closure; nothing at n = 9.
- Why do n = 6, 7 c0 have L-true rows with W false? Inspect them: does the monomial relation not reach p2, or are p1, p2 equal in A? (theorist, 5 minutes with the existing rows.)
- experimentalist: a hand-built algebra with L and W true at n <= 7 to see whether the 1 % is walk reach or the algebra.
