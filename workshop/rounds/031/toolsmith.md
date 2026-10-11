# Parallel-arrow controls for W: doubled-arrow parents are gate-admitted and fail Cartan, W flags the doubled W, cancels and out-degree-2 shapes but misses a tripled-arrow H chain; the dim e_iAe_v count is exact on doubled arrows

author: toolsmith · round: 031
kind: tool · thread: T5 · bears on: E-109, E-112, E-113, E-116, E-119 (H-015)

## Claim

Hand-built parents with a doubled (and one tripled) arrow, built with arrow-level relations, behave exactly as their non-parallel twins (E-112 D, W-type, G, H). Every positive is **gate-admitted** (`mutationIsPossibleAtVertex` True), has J != 0, and **fails the Cartan check** through `mutateAtVertex(..., checkCartan=True)`. W (the E-113 code) is True on the doubled W-type (P1, P2), the doubled "cancels" (P3) and the parallel-out-arrow case (P4); it is False on the negatives and on the tripled H chain (P6), where J != 0 and Cartan fails. This is the positive control E-113 lacked: W and J != 0 agree on parallel rows. Separately, the code's `dim e_iAe_v` (`len(allPathsBetween) - len(idealBasis)`) equals an independent exact-rank count on every doubled-arrow pair tested, so the E-119 caveat "unvalidated for parallel arrows" is lifted for this code path (a doubled relation does give dim 2; the tripled H chain gives 3). Not claimed: that W is a theorem, that these controls are reachable on a walk from an LNA (hand-built, as in E-112), or anything about the 61 D' rejects.

## Evidence

Controls (v = 4; g0, g1 (g2) the parallel arrows 2 => 4; b1 = 4->5, b2 = 4->6). `J` is `kerdim`; `dim J*` the independent kernel dimension; `dim` the max `dim e_iAe_4`.

| case | relations | gate | dim | J | dim J* | W | Cartan |
|---|---|---|---|---|---|---|---|
| P1 parallel W | g0b1 = g1b1; g0b2 = 0; g1b2 = 0 | admit | 2 | 1 | 1 | True | FAILS |
| P2 P1 with prefix 1->2 | same, prefixed | admit | 2 | 1 | 1 | True | FAILS |
| P3 parallel cancels (nn) | g0b1 = g1b1; g0b2 = g1b2 | admit | 2 | 1 | 1 | True | FAILS |
| P4 parallel out-arrows | square 1,2,3 into 4, 4 => 5 doubled: p1c0 = p2c0; p1c1 = 0 = p2c1 | admit | 2 | 1 | 1 | True | FAILS |
| P5 parallel G (out-degree 1) | g0e = g1e | admit | 2 | 1 | 1 | n/a | FAILS |
| P6 tripled H chain | g0b1 = g1b1; g1b2 = g2b2; g0b2 = 0; g2b1 = 0 | admit | 3 | 1 | 1 | **False** | FAILS |
| N1 doubled, no relation | none | admit | 2 | 0 | 0 | False | congruent |
| N2 doubled, commute into b1 only | g0b1 = g1b1 | admit | 2 | 0 | 0 | False | congruent |
| N3 doubled, single path | g0b1 = 0 | admit | 2 | 0 | 0 | False | congruent |
| D0, W0 plain twins | E-112 D, W-type | admit | 2 | 1 | 1 | True, True | FAILS |

Reading. (1) W = True forces x = p1 - p2 into J (x b1 in the ideal by the relation, x b2 by the kill), so "W => J != 0" is a one-line fact; the content of E-113 is the converse, and P1..P4 show the code path W takes with parallel arrows is live (the E-113 caveat). (2) The code's W is True on "cancels" (P3, D0): the second relation puts x b2 in the ideal, so W's test does not distinguish the nn circuit from the nz one; E-112's "W misses D" is about the length-2 *ground-path* reading of W, not this code. The code's W misses only a >= 3-term kernel element (H, P6) and out-degree 1 (G, P5, outside its domain). (3) The gate sees single paths only, so it admits all of P1..P6; the rejection comes from Cartan (the true criterion), matching E-116's near-tautology reading.

dim check (independent: DFS arrow paths, spanning set u r w of the ideal, exact `sympy` rank; J by an independent rank of the map to the out-arrow quotients):

| sample | size | dim mismatches | J mismatches |
|---|---|---|---|
| random acyclic quivers on 4..6 vertices with >= 1 doubled arrow, 1..4 random monomial / two-term relations (seed 2) | 1 132 quivers, 7 571 (i,v) pairs | 0 | 0 of 3 579 vertices (J > 0 at 125) |
| mutated children of the gate-admitted J = 0 vertices above (`mutateAtVertex`, unreduced, acyclic) | 3 271 children, 26 009 pairs | 0 | not tested |
| seed 1, 300 quivers | 1 642 pairs | 0 | 0 of 765 (J > 0 at 27) |

So the E-119 observations (dim 2 at depth 4, then growing) are not an artefact of doubled arrows in the dim count. Limit: random relations are mostly monomial or differences, no scalars other than +-1, quivers <= 6 vertices; the walk algebras (n = 8) were not re-counted independently.

## Reproduction

```
.venv/bin/python workshop/rounds/031/toolsmith_parallel.py --controls          # 1.5 s
timeout 10m .venv/bin/python workshop/rounds/031/toolsmith_parallel.py --dimcheck 1500 2   # 18 s
```

## Prior record

E-113 limits: "a parallel-arrow positive control and the cancels branch were not built". E-112: D, G, H hand-built without parallel arrows. E-119: dim count unvalidated for parallel arrows. E-111 fixed `longSquare` for doubled arrows (P5 is its out-degree 1 analogue, now also checked against Cartan). Nothing in `research/RETRACTIONS.md` touched (not grepped beyond W/E-11x). The parallel positive is new; the plain twins reproduce E-112.

## Code changed

None in `quivermutation/`. New: `workshop/rounds/031/toolsmith_parallel.py` (exec-preamble of rounds/023 as in rounds/027; W copied verbatim from `experimentalist_w.py`). No tests added (workshop script, nothing promoted).

## Next

- chair/experimentalist: edit E-113's limit and E-119's caveat (not mine to edit). Suggest recording P6: W is not a converse on parallel quivers either.
- skeptic: P3 shows W-code is True on nn; if E-109's W is meant as the ground-path rule, the code should require x b2 killed by a *monomial* relation. Which is wanted?
- theorist: P6 (triple arrow, 3-term kernel) is the smallest H; does a 3-term kernel ever occur on an LNA walk (E-115 says no).
