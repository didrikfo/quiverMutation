# The E-174 End(T) comparison rejects 280 of 280 equal-dimension, non-isomorphic parallel-arrow children (n = 7, class 1), so its 25/25 "iso" is not vacuous on the relations -- but it still cannot see J, and takes dim as an unchecked input

author: skeptic · round: 058 · kind: result
thread: T10 · bears on: E-174, E-167, E-168 (no H-/F- status)
scope: n = 7, class 1 only (E-151 walk rebuilt at 20 000 expansions); the 7 failing steps whose child has exactly one parallel pair (7 of the 8 parallel-arrow ones; step 12 has two pairs, not done) and exactly 3 line-type relations; one family of wrong algebras (re-chosen killed line); label-preserving iso only; class 2 and n = 10 not done.

## Claim

For each of those 7 steps, I replaced the child's three "line-type" relations (x then the parallel pair P = (b0,b1) kills one combination of x b0, x b1) by every choice of killed line among [1:0],[0:1],[1:1],[1:-1] (4^3 = 64 variants per step, 448 in all). Every variant has the SAME Cartan matrix, SAME quiver, and the same dim K Q/I' as the true child (dimension recomputed here by linear algebra, not read from the test). 280 of the 448 are not isomorphic to the true child (the three killed points of P^1 have a different coincidence pattern; a Moebius map preserves it), 168 are. The E-174 comparison (`symcheck2`) answers "NO" on all 280 and "iso" on all 168: 448/448 correct. So at n = 7 class 1 a wrong parallel-pair algebra with equal dimension vector is rejected; the E-174 test has power against relation errors of this kind.

It does NOT claim: power against a wrong algebra of different dimension (see Evidence 2: the test never compares dim K Q/I' with dim End(T); it takes the child's dims as input and accepts a surjection), against non-label-preserving or non-parallel wrong algebras, at n = 8/10, or in class 2. Refuted by any equal-dim non-iso variant accepted, or an iso variant rejected.

What the 25/25 + 370/370 + 35/35 therefore show: given that the child's relation set `crels` and dimension are right (inputs), End(T) as an algebra is that child, label-preserving, including the relation structure at parallel pairs. They do not show anything about J: End(T) is iso to the child at J != 0 steps too (E-174), so the comparison is blind to J by construction (Hom(T,T[-1]) != 0 is what separates). The 25 failing steps' derived-equivalence status still rests on Hom(T,T[-1]) = 0 at the J = 0 path edges and E-168 generation, not on this comparison.

## Evidence

1. Table, equal-dimension variants only (the 448 are all equal-dim, `quotient_dims` agrees with the child's Cartan on every pair):

| truth (independent) | test = iso | test = NO |
|---|---|---|
| iso to true child | 168 | 0 |
| not iso | 0 | 280 |

Steps 0, 1, 2, 3, 8, 10, 11 (c1 numbering of E-174); true killed points per step printed by the script. With three points on P^1 the true algebra has all three lines distinct; the 280 wrong ones have a coincidence among the three (or reorder with a different coincidence pattern). Note a limit of the family: three points have no moduli, so the control separates "pattern" only, not a continuous parameter (a 4th line would; none occurs here).

2. Dimension input. Dropping one relation (47 variants over the parallel-arrow steps): dim K Q/I' > dim End(T) in 47 of 47, `symcheck2` says "iso" in 47 of 47. The proof of "onto => iso" in the E-174 docstring needs equal finite dimensions; the code checks only arrow counts per pair and, in `compare2`, Cartan equality of the child with End(T) (true child dims). So the soundness of "iso" depends on `crels` being a complete relation set for the child -- exactly the input E-174's Scope names. This control could not move that: I perturbed `crels`, not the child algebra as produced by `quiverMutationAtVertex`/`procedure.relationsFrom`.

3. Consistency with record: E-174 reports "binomial -> first term rejected 9/9, doubled coefficient absorbed 9/9". My family is the same kind of perturbation made exhaustive and with independent dimension and truth; the one new fact is that the rejected perturbations were equal-dim and genuinely non-isomorphic (E-174 did not check that), and that equal-dim non-iso is rejected 280/280.

## Reproduction

```
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100   # twice; ~10 min total
CLS=1 timeout 10m .venv/bin/python workshop/rounds/058/skeptic_wrongalg.py     # ~1 min; output tail in skeptic_wrongalg_out.txt
```

## Prior record

E-174 Scope lists "an equal-dims wrong-algebra control" as not covered; E-174 power: 9/9 and 9/9 (unchecked dims). Nothing in RETRACTIONS bears on this. New: the control, its 448-row table, and the dimension-input observation (the 47 drops).

## Code changed

None to library; new `workshop/rounds/058/skeptic_wrongalg.py` (exec's `toolsmith_endt2.py`). No tests touched.

## Next

- toolsmith: have `compare2` compute dim K Q/I' from `crels` and assert it equals dim End(T), closing the dimension input (cheap: `quotient_dims` in my script).
- skeptic/toolsmith: the same table for step 12 (two parallel pairs, where a 4th-line modulus could appear) and class 2; a wrong algebra from the actual search (a non-key-keeping neighbour) rather than perturbed relations.
- The remaining T10 premise is unchanged: J = 0 along the paths (Hom(T,T[-1]) = 0) and the completeness of `crels`.
