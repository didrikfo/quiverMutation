# Peeling works for LNAs with several big relations once chain links are relations of any length, ending at x1 on the left and starting at x_{m-1} on the right; fits n = 8..10 exhaustively and 8 of 8 new predictions at n = 10, 11

author: maverick · round: 023 · kind: result
thread: T6 · bears on: E-103, E-101, H-017

## Claim

Take LNAs with two or more relations of >= 3 arrows ("big"). Write a relation as (s, e), e = s + length (zero path x0 = s .. xm = e).
**D1** (E-103: depth 1 iff some big relation is not blocked) holds on all of them, n = 8 (236 LNAs) and n = 9 (981): nothing new, it was
already 0 mismatches over every LNA at n = 6..10. The **peeling formula as recorded** (links = 2-arrow relations only) fails exactly at
the LNAs whose blocker is itself a big relation: 1 at n = 8 (`230302`), 6 at n = 9 (`0230302 2230302 2303020 2303022 2304002 2400302`),
24 at n = 10; each has true depth 2, formula 1. In each, a relation of 3+ arrows blocks one end of a big relation (right end of s = 2 in `230302`:
the relation starting at x_{m-1} = 4 has 3 arrows, so the old chain count b = 0, though D1 counts it as blocking).

**Mirror-chain formula (a fit, not a proof).** For a big relation (s, e): left chain = longest sequence of relations R1, R2, ... of any
length with e(R1) = s + 1 and e(R_{k+1}) = s(R_k) + 1; right chain = longest sequence with s(R1) = e - 1 and s(R_{k+1}) = e(R_k) - 1.
depth = 1 + min over big relations of min(left, right). For 2-arrow links this is E-103's formula (both sides reduce to consecutive
2-relations), so single-big-relation LNAs are unchanged. It is the formula D1 itself suggests: D1's "blocked" tests for any relation ending at x1 / starting at x_{m-1}, not a 2-arrow one.

**Cord depth when the old formula fails:** 2 for every one of the 7 LNAs at n = 8, 9 (and 30 at n = 10 with depth 2); the one case with two big relations needing 3 is `22303022` at n = 10, and 7 at n = 11 (below).

Not claimed: a derivation of the chain step; any statement for n >= 12; that the old formula is "wrong" for single-big LNAs (it is the special case).

## Evidence

All numbers use the recorded depths (`rounds/022/theorist_blocked_depths.txt`, every LNA with depth >= 2 at n = 8..10; every LNA absent from it has depth 1 by D1).

| n | class | LNAs | old formula mismatches | mirror-chain mismatches |
|---|---|---|---|---|
| 8 | one big | 129 | 0 | 0 |
| 8 | two or more big | 236 | 1 | 0 |
| 9 | one big | 321 | 0 | 0 |
| 9 | two or more big | 981 | 6 | 0 |
| 10 | one big | 769 | 0 | 0 |
| 10 | two or more big | 3837 | 24 | 0 |

(E-103's "31 mismatches" = 1 + 6 + 24 over n = 8..10; matches.) The fit was made on these data, so the n = 10 row is not independent of the n = 8, 9 rows only in that
the rule was read off the n = 8 failure; the real test is out of sample. The chain formula predicts depth 3 for two-big LNAs where the old one says 1:
`22303022` (n = 10) and seven at n = 11 (`022303022 222303022 223023022 223030220 223030222 223040022 224003022`; n = 12 not enumerated).
Search with `theorist_cordcrit.py` (the E-103 visitor): L = 2 finds no cycle member in all 8, L = 3 finds one in all 8, so depth is exactly 3 (8 of 8,
nothing one short; seconds to 33 s per batch). Predicted-depth distribution of two-big LNAs (chain): n = 8 {1: 235, 2: 1}, 9 {1: 974, 2: 7},
10 {1: 3806, 2: 30, 3: 1}, 11 {1: 14379, 2: 105, 3: 7}.

Weakness: all failing two-big cases at n <= 10 have depth 2, so the n = 8, 9 two-big data cannot tell the mirror chain from "any blocked LNA has depth >= 2";
the discriminating evidence is the 8 depth-3 predictions only, all of one shape (two 3-relations in contact, `303`, flanked by 2-chains). Two big relations
not in contact never matter in these data (the min is attained by the unblocked one). Depth data are minimum cycle-member depths of the repo's mutation search, and "cycle member" is not the GLOSSARY cord (E-103).

## Reproduction

```
.venv/bin/python workshop/rounds/023/maverick_twobig.py 8        # old formula vs data, two-big LNAs (seconds); also 9
.venv/bin/python workshop/rounds/023/maverick_variants.py        # four variants of the old formula, n = 8..10, all fail the same LNAs
.venv/bin/python workshop/rounds/023/maverick_chain.py           # mirror chain: 0 mismatches n = 8..10 (about 1 min)
.venv/bin/python workshop/rounds/023/maverick_predict2.py        # predicted depth >= 3 two-big LNAs at n = 10, 11
NAMES=22303022 timeout 10m .venv/bin/python workshop/rounds/022/theorist_cordcrit.py 10 3      # depth 3 (7 s; L = 2: none)
NAMES=022303022,222303022,223023022,223030220,223030222,223040022,224003022 timeout 10m .venv/bin/python workshop/rounds/022/theorist_cordcrit.py 11 3   # outputs maverick_n11_L2.txt, _L3.txt (28 s, 33 s)
```

## Prior record

E-103 states the formula, the 31 mismatches and "no formula" for two big relations; this supplies one. E-101 (n = 8) is unchanged. No other entry
mentions it (grep `230302`, "mirror chain").

## Code changed

None to existing files. New: `maverick_twobig.py`, `maverick_variants.py`, `maverick_chain.py`, `maverick_predict2.py` and the two n = 11 outputs, all in `workshop/rounds/023/`.

## Next

- Skeptic: pick the n = 12 two-big LNAs the chain formula puts at depth 4 (and long-chain ones with a big relation as link at depth >= 3 that are not of the `303` shape), run them. One miss kills the any-length link rule.
- Theorist: derive why a blocker of any length on x1 / x_{m-1} costs one step and the next link is at start+1 / end-1; why `303` costs the same as `202`.
- Toolsmith: a fast depth-1 cycle test (D1 on the sequence) would let n = 12 be enumerated without the search.
