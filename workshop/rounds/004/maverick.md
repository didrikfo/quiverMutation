# H-017 survives to depth 6 at n = 9, but the Euler form / Coxeter polynomial cannot predict it; what the Euler form does give is a signature test for "outside every quipu class"

author: maverick · round: 004 · kind: result (H-017 tested; reframing mostly negative)
thread: T6 · bears on: H-017, F-034, F-045, F-047, F-048

Speculation level: **tested on small cases** (n = 9 for H-017; n = 8..11 for the signature test).

## Claim

1. **H-017 holds on everything I could reach, and not because of anything the invariants force.** Over the 9 LNAs at n = 9 outside a quipu class, the quipu-with-relations algebras *proved* in their classes (mutation walk to depth 4 for all nine, 5 for all nine, 6 for `3033030` and `4444400`) have (cords, relations) exactly the ten pairs of the record, and none has relations <= cords (0 of 45..97 per LNA at depth 6). Not refuted, not proved. Depth 6 does not move the per-LNA minimum of `3033030` (2); it does move it for `4444400` (2 at depth 5, 1 at depth 6), which is in the 8-member class whose minimum was already 1, so the "class invariant" reading in H-017 is fine but per-LNA minima are depth-dependent.
2. **The Coxeter polynomial does not see (cords, relations).** Every tree algebra with monomial relations has `c_{n-1} = 1` (trace of Phi is -1: arrows contribute `n-1`, a zero relation contributes 0 to `sum C_ij (C^-1)_ij`), and `c_{n-2}` is not a function of (cords, relations) (n = 6, table below). At n = 9 the polynomial matches 3677 quipu-with-relations algebras (`families.py quipus 9`), of which 1465+190 fail a Smith-form test and the rest include many below the diagonal. So the polynomial-only candidates below the diagonal (the "sharp prediction" of H-017) are not excluded by the polynomial, by the Smith normal form of `C + C^T`, or by the full F-047 profile (the two filters keep identical survivors).
3. **What the Euler form does say, and the plain census does not.** The signature of `G = C + C^T` (congruent to `C^-1 + C^-T`, so a derived invariant) separates the classes outside every quipu class at every length checked: for all LNAs of length 8, 9, 10, 11, `pos(G) <= n - 2` **iff** the LNA lies in no quipu class (11: 2647/2647 non-quipu incl. 16 unplaced have `pos = 9`; all 14 149 quipu-class LNAs have `pos >= 10`; 9 and 10 likewise). One direction is a proof by computation for `n <= 12`: a hereditary quipu has `G = 2I - A`, and no quipu of order <= 12 has two adjacency eigenvalues >= 2 (127 quipus at 12, checked), so its class has `pos >= n-1`. It fails first at n = 13 (`P^(1,0,3,0,1)_(1,1,1,1)`, two disjoint D~ pieces; 1, 11, 74 quipus at 13, 14, 15). Consequence for relations: one relation (gldim <= 2) is a rank-2 perturbation of `G` with one positive and one negative eigenvalue, so it lowers `pos` by at most 1; hence such a quipu needs at least `2 - b` relations, b in {0,1} the number of adjacency eigenvalues >= 2. That is "relations >= 1 or 2", independent of cords. It cannot give "relations > cords" (which is 4 relations at 3 cords), and I did not try to make it.

It does **not** claim: that relations > cords is a theorem, that below-diagonal candidates are outside the class (the direct search below is weak), or that the signature criterion is a characterisation beyond n = 12.

## Evidence

Census of polynomial-match candidates at n = 9 by (cords, relations), counting orientations and ideals as `families.py` does, after the Smith-form filter (survivors only). rels - cords < 0 in bold.

| class polynomial | survivors | below diagonal | cells (cords,rels): count |
|---|---|---|---|
| `3033030` | 596 | **(3,1) 36, (3,2) 24** | (1,3)10 (1,4)12 (1,5)12 (2,2)33 (2,3)65 (2,4)106 (2,5)32 (2,6)12 (3,3)142 (3,4)28 (3,5)84 |
| other 8 LNAs (tubular) | 1426 | **(2,1) 24, (3,2) 76** | (1,2)20 (1,3)84 (1,4)94 (1,5)46 (1,6)4 (2,2)137 (2,3)258 (2,4)243 (2,5)52 (2,6)12 (3,3)160 (3,4)156 (3,5)60 |

Direct test of the sharp prediction: 16 below-diagonal survivors (4 per cell, cells (3,1), (3,2) x2 shapes, (2,1)), mutation search to depth 4 from each (both orientations of the search): reached **no** LNA of the class. Weak: no positive control (several diag >= 1 candidates I tried at depth 5 also reached nothing), so this neither supports nor refutes it.

Walk data (union of pairs over all LNAs; depth 4 reproduces the record exactly): depth 5 adds nothing new to the ten pairs; depth 6 for `3033030`: (1,3)(1,4)(1,5)(2,4)(2,5)(2,6)(3,5); for `4444400`: (1,2)(1,3)(1,4)(1,5)(2,4)(2,5)(2,6)(3,4)(3,5).

Second coefficient, all quipu algebras with relations at n = 6 (c_{n-1} = 1 in all): (cords,rels) = (1,1) gives c_{n-2} in {1, 0, -1}; (2,2) in {-1, 0, 1}; (2,3) in {0, -2}; (2,4) in {-1, -5}. Not a function.

Signature `(pos, neg, zero)` of `C + C^T` by status, all LNAs:

| n | quipu class | not quipu / unplaced |
|---|---|---|
| 8 | (7,0,1) 77, (7,1,0) 90, (8,0,0) 262 | none |
| 9 | (8,0,1) 415, (8,1,0) 733, (9,0,0) 273 | (7,0,2) 9 |
| 10 | (9,0,1) 104, (9,1,0) 3919, (10,0,0) 577 | (8,0,2) 129, (8,1,1) 130, (8,2,0) 1; unplaced (8,2,0) 2 |
| 11 | (10,0,1) 272, (10,1,0) 12660, (11,0,0) 1217 | (9,0,2) 25, (9,1,1) 1848, (9,2,0) 758; unplaced (9,1,1) 16 |

## Reproduction

All from the repository root, `.venv/bin/python`, prefix `timeout 10m`.

```
workshop/rounds/004/maverick_census.py 9                # candidate cells by (cords, rels), 7 s
workshop/rounds/004/maverick_reached.py 9 4             # proved members' pairs, all 9 LNAs, 107 s
workshop/rounds/004/maverick_reached.py 9 5 4550400 5504030 5040330 5050030   # ~8 min; the other five LNAs alone, ~2 to 5 min each
workshop/rounds/004/maverick_reached.py 9 6 3033030     # under 10 min; same for 4444400
workshop/rounds/004/maverick_euler.py 9                 # signature + Smith form of C+C^T vs candidates, ~1 min
workshop/rounds/004/maverick_profile.py 9               # F-047 profile filter (same survivors)
workshop/rounds/004/maverick_verify.py 9 4 4 -1         # 16 below-diagonal candidates searched, 195 s
workshop/rounds/004/maverick_coeff.py 6                 # c_{n-2} not a function of (cords, rels)
workshop/rounds/004/maverick_corank.py 9 (or 8, 10, 11) # signature by status; n = 11 takes 4 min
workshop/rounds/004/maverick_quipu_pos.py 15            # quipus with two adjacency eigenvalues >= 2
```

## Prior record

H-017 and its census are `EXPERIMENTS.md` (E-030 area, section 5, "What the quipu members look like") and F-034. Euler form ideas already recorded: F-045 (derived-tame iff Euler form PSD; corank + Dynkin type classify tame LNAs; 3033030 is corank 2, type D_7), F-048 (periodic Coxeter + indefinite form certifies not PH). H-017 itself points at the Euler form but records no test of it. I found no record of the signature criterion "`pos(G) <= n-2` iff outside all quipu classes" (grep signature/corank/positive/Euler over `research/`, nothing beyond F-045/F-048), nor of the n = 13 breakdown; treat as new but modest: it extends F-045/F-048's "Dynkin and Euclidean have no nonpositive directions beyond one" to "quipus have at most one". For n = 9 it is F-045 restated: the 9 outside LNAs are precisely the corank-2 ones.

## Answer to "does the reframing say anything the plain version does not"

Two things, both small. (a) It explains *why* a quipu class needs relations in a member with the right quiver shape (tree-lattice signature has `pos >= n-1` for n <= 12), a bound of 1 or 2 relations, independent of cords. (b) It gives a cheap exact test for "outside every quipu class" for n <= 12, without any walk. It does **not** predict cords or relations; the polynomial and Cartan-congruence data are blind to them; the Coxeter polynomial trace is the same for every tree algebra. So the reframing does not explain or test H-017's "more relations than cords".

## Code changed

None in the library. Scripts only, all in `workshop/rounds/004/` (`maverick_*.py`); no tests needed.

## Next

- Overnight proposal: `maverick_reached.py 9 7` for `3033030` and `4444400`, and `maverick_reached.py 10 4` over the 262 outside LNAs (H-017 at n = 10), each shard by class name as the script allows; below-diagonal pairs are the only outcome that matters.
- A positive control for the direct search: take a certified member (from `reachedQuipuAlgebras`) and check `verify`-style search from it finds the LNA at its path depth, else the 16-candidate negative means nothing (toolsmith).
- Theorist: prove or refute `relations >= cords + 1` in the gldim-2 case by an eigenvalue argument that uses the tree (each cord is a branch vertex, adding an adjacency eigenvalue of the tree); my rank-2 count only gives 1 or 2.
- Theorist/experimentalist: check the signature criterion on LNAs at n = 12 (58 786 rows, ~15 min) and at n = 13, where a quipu class with `pos = n-2` should first appear.
