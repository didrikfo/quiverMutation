# The 25 loose-shape W-false rows at n = 6, 7 are all half-W (one term killed by the monomial, the other not); the gate and W do not differ from n = 8, what differs is that n = 6, 7 walks contain no both-terms-die row

author: skeptic · round: 031 · kind: result
thread: T5 · bears on: E-109, E-113, E-115, E-116

## Claim

All 25 loose-shape (L) rows at n = 6, 7 class 0 (5 + 20, rows of E-116) fail W in the same way: for R = p1*b1 = p2*b1 and the monomial M ending in b2 at v,
exactly one of p1*b2, p2*b2 is zero in A (the term that has M as a suffix; the other term is not killed). Hence x = p1 - p2 has x*b2 = the surviving term != 0, x is not annihilated by b2,
and J = 0 (K false in all 25). So the answer to "W or the gate?": neither differs. W is the same rule, the gate (single-path test) is the same, and J = 0 agrees with W-false row by row.
What differs is the population: the 25 rows are all "arrow versus path" squares (term lengths (2,1) or (1,2): one side is a single arrow into v, which M kills outright), whereas the 61 rejects have both terms of length >= 2 and both killed.
It is NOT a size effect of the L shape: the 23 L-true accepts at n = 8 c0 have the identical half-W pattern, and include lengths (2,2) and (2,3). Not claimed: why n = 6, 7 walks never produce a both-die row (walk reach versus algebra), nor that none exists at n = 6, 7 (prefix walks, capped).

## Evidence

Classification per (relation R, monomial M) pair, same walk and caps as E-116 (n = 6: 3 000 algebras, n = 7: 2 522, n = 8 c0: 5 145). Pattern = (p1*b2 in I, p2*b2 in I). "Half" = exactly one in I; "Both" = both in I. In every row the W-clause (p1 - p2)*b2 in I equals "Both", and K equals "Both".

| n, c0 | L pairs | Half (W false, K false) | Both (W true, K true) | lengths of Half rows (p1,p2) | M length |
|---|---|---|---|---|---|
| 6 | 5 | 5 (3 with p2 the short term, 2 with p1) | 0 | (2,1) x3, (1,2) x2 | 2 |
| 7 | 20 | 20 (8 + 12) | 0 | (2,1) x8, (1,2) x12 | 2 |
| 8 | 84 pairs (reported as 84 rows in E-116; here 84 = 23 + 61) | 23 | 61 | (2,1)/(1,2) x14, (2,2) x5, (2,3)/(3,2) x4 | 2 |

Two sub-patterns in the 25: (a) M = the arrow-pair (p_short, b2), p_short a single arrow into v, long side p_long = c d with d into v and (d, b2) nonzero: 100% of the 25.
Example n = 6: R = (3,1)(1,6)(6,2) = (3,6)(6,2), M = (3,6)(6,4): the direct arrow 3->6 dies into 6->4, the route through 1 does not. For the rejects, both sides pass through different mid vertices and
each side meets its own zero relation (n = 8 example: (4,2)(2,8)(8,3) = (4,6)(6,8)(8,3) with a zero relation on each side into (8,1)), i.e. two independent monomial-type facts, not one.
M-suffix flags: in all Half rows exactly one term has M as a suffix; in Both rows the second term dies by a different relation (not captured by my one-M flag, so this half of the mechanism is "dies in A", not "dies by M").

Consequence for the gate: the two reasons are the same fact. The gate refuses v iff some nonzero single path q into v has q*beta = 0 for all out-arrows; in a Half row neither p1 nor p2 qualifies (b2 does not kill p1 or p2 simultaneously), and x is not killed by b2, so the gate admits and J = 0. Half rows are genuine accepts, not a gate failure.

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/031/skeptic_loose.py 6 0 3000   # ~100 s
timeout 10m .venv/bin/python -u workshop/rounds/031/skeptic_loose.py 7 0 2522   # ~110 s
timeout 10m .venv/bin/python -u workshop/rounds/031/skeptic_loose.py 8 0 5145   # ~170 s solo
```
Outputs `workshop/rounds/031/skeptic_loose_n{6,7,8c0}.txt`, with one example relation pair per pattern.

## Prior record

E-116 left this open; E-115 says the 22 rows are "half-W (20) or two loose pendants (6)" -- a different set (J nonzero rows) with a similar word. E-109/E-113: W has 0 mismatches; consistent here (K = Both in all 109 pairs). New: the 25 n = 6, 7 rows are exactly the n = 8 accepts' pattern, not a separate phenomenon. Nothing in RETRACTIONS touched.

## Code changed

None in the library. New `workshop/rounds/031/skeptic_loose.py`; no tests.

## Next

- Weak spots: counts are per (R, M) pair, E-116's 5/20/84 happen to equal them; my sub-claim (a) "100% of 25 have a single-arrow term" is read from the table, the example per pattern only is printed; the Both-row second-kill mechanism is not classified.
- experimentalist: hand-build an n = 6 both-die L algebra (two mid vertices, a zero relation on each side into b2) and test gate admission and reachability by the walk; if the algebra is a valid parent but never reached, the n <= 7 absence is walk reach.
- theorist: is there a reason the walk from LNA keys never yields a (>=2, >=2) both-die square below n = 8 (vertex count 6 suffices in the square itself)?
