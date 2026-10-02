# The skeptic's 22 rows are not length-2 ground paths (half-W or two loose pendants, no circuit); on walks every Gamma_i component has at most 2 edges, and the nn / long-circuit question reduces to dim e_iAe_v >= 2 / >= 3, which is rare, but is NOT proved

author: theorist · round: 029 · kind: result + negative (no proof of the obstruction) · thread: T5 · bears on: E-105, E-107, E-108, E-110, E-111

## Claim

1. **The 22 rows.** At n = 8 class 0 (150 s walk, 8 030 algebras, 3 716 out-degree 2 rows) there are 28 (row, i) pairs with J_b1 != 0, J_b2 != 0 and J = J_b1 ∩ J_b2 = 0 (the skeptic's 22 are the row-level version; cap-dependent). For 26 of the 28 the circuit graph Gamma_i has a single nontrivial component of exactly two edges, and it is **not** a ground path (it has no circuit: J = 0 forces that):
   - **20 "half-W"**: edges (a, c) and (a, g): classes p1, p2 with p1 b1 = p2 b1 != 0 (the commutativity into b1 that W needs), and p2 b2 = 0 but p1 b2 != 0. The path starts at ground and ends at the nonground vertex [p1 b2]. So J_b1 = <p1 - p2> and J_b2 = <p2>, which meet in 0. This is W with one of the two kills missing.
   - **6 "two loose pendants"**: edges (a, g), (g, c) joined only through ground: p1 b2 = 0 and p2 b1 = 0, no shared product. J_b2 = <p1>, J_b1 = <p2>, two different single paths; each arrow kills a different path, no relation links them.
   - 2 are outside the lemma (a relation with a repeated term, i.e. a coefficient 2: `5 4 3 + 5 4 3 + 5 7 3`; the classes are not single products). Unresolved, not claimed.
   So the answer to "single length-2 ground path?" is **no**: the 22 are the near-misses of W and of the length-2 ground path, which is what the skeptic's control suggested (the 42 have both kills, the 22 have one, or one each on different paths). The same census at n = 8 c1 (13 621 algebras, 6 283 rows): 0 such pairs.
2. **Every J != 0 pair is the length-2 ground path** (46 at c0, 10 at c1; shape (edges, nonground vertices, ground endpoints) = (2, 1, 2)), re-confirming E-110 with the graph computed directly (class independence and single-class products checked: 0 failures except the 4 coefficient-2 pairs above).
3. **Component size.** Over all (row, i) with a path i -> v (c0: 6 776 pairs, c1: 8 990), every component of Gamma_i with at least 2 edges has exactly 2 edges (shapes seen: (2,1,2) 46, (2,2,1) 20, (2,2,2) 7, (2,3,0) 11 at c0; (2,1,2) 10 at c1). A component with >= 2 edges forces J_b1 or J_b2 != 0 at that i, so the J, 22-type and one-sided categories list all of them (84 pairs at c0, 10 at c1); the other pairs have only single-edge components. **No nn 2-cycle (shape (2,2,0)), no circuit of length >= 3, no component of >= 3 edges (ground counted as one vertex) occurs.**
4. **Reduction (proved, easy).** A circuit with k edges uses k distinct classes of e_iAe_v, linearly independent (checked). So a circuit of length >= 3 needs dim e_iAe_v >= 3, an nn 2-cycle needs dim >= 2 with p1 b = p2 b != 0 for both b. Observed dim e_iAe_v over admitted out-degree 2 vertices: c0 {0: 1020, 1: 5492, 2: 261, 3: 3}; c1 {0: 1448, 1: 7480, 2: 62}. Only 3 pairs ever have dim 3 and none carries a long circuit; none of the 323 dim-2 pairs is nn.
5. **Obstruction by key for the hand cases (checked).** D (nn), the W-type control, G and H have Coxeter keys (1,2,-1,-4,-1,2,1), (1,1,-2,-4,-2,1,1) (both W-type and G), (1,3,7,13,13,7,3,1) which are not keys of any LNA or dual LNA with the same number of vertices (6, 6, 6, 7). So none of the four is derived equivalent to an LNA and cannot occur on a walk. This disposes of those four examples only; it is **not** an obstruction for a hypothetical other nn or long-circuit algebra, which would have to be built to land on an LNA key.

What is not claimed: that no nn 2-cycle or circuit >= 3 can occur on LNA-derived algebras (no proof; the evidence is a capped BFS of two classes at n = 8); anything at n = 9, classes 2-3, out-degree >= 3, scalars != 1 (4 pairs out of lemma), or parallel arrows (their products are not separated in this script's class basis; same code path as `scholar_pairtest`).

## Evidence

Script `theorist_circuit.py` rebuilds the round-023 walk (cap 150 s, so counts differ from the skeptic's 200 s run: 7 848 algebras there), for every out-degree 2 gate-admitted row and every source i: J, J_b1, J_b2 by rank, the classes of e_iAe_v (normal forms up to scalar), their products with b1, b2 (zero = ground), then connected components with ground as one vertex.

| walk | algebras | rows | (row,i) J != 0 | 22-type | one-sided | max dim e_iAe_v | components with >= 3 edges |
|---|---|---|---|---|---|---|---|
| n=8 c0 | 8 030 | 3 716 | 46 | 28 (26 analysable) | 2 096 (2 not analysable) | 3 | 0 |
| n=8 c1 | 13 621 | 6 283 | 10 | 0 | 2 781 | 2 | 0 |

22-type shapes at c0: (2,2,1) 20, (2,2,2) 6. Weakest steps: (a) "ground counted as one vertex" merges components that a circuit analysis keeps apart; I used it only to bound component size, where it can only make components larger, so the bound is safe; (b) the class basis check (independence, single-class products) is by normal form, trusted for scalar-1 relations only; (c) the 4 coefficient-2 pairs.

Why the 22 exist, in one line: Gamma_i has the circuit shape of W on the b1 side (shared product) but the b2 side only grounds one of the two paths. Whether the Nakayama relation shape forces "at most one of the two paths is killed by b2" (hence at most half-W) is the open mechanism; the key obstruction (item 5) is the only proof-level statement and it is example-specific.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/029/theorist_circuit.py 8 150 0    # about 3 min, output theorist_circuit_n8c0.txt
timeout 10m .venv/bin/python workshop/rounds/029/theorist_circuit.py 8 150 1    # about 3 min, output theorist_circuit_n8c1.txt
timeout 10m .venv/bin/python workshop/rounds/029/theorist_keys.py               # 30 s, the Coxeter keys of D, W-type, G, H
```
Counts move with load (wall-clock cap); the shapes should not.

## Prior record

E-110 (circuit lemma, D/G/H, "on walks always a length-2 ground path", why no nn or longer is open), E-107/E-108 (W; the 42), the skeptic's 22 (round 026 `skeptic_x_out.txt`, "not explained"). Not recorded before: the shape of Gamma_i for the 22 (half-W and loose pendants), the component-size bound 2, the dim e_iAe_v tally, and the key obstruction for D, G, H (E-103 only says two examples have Coxeter polynomial in no class at n = 6; I did not check which). Nothing in RETRACTIONS bears on it (grep "circuit").

## Code changed

New: `workshop/rounds/029/theorist_circuit.py`, `theorist_keys.py`, outputs `theorist_circuit_n8c{0,1}.txt`. No library change, no tests run.

## Next

- Theorist: dim e_iAe_v = dim Hom_D(T_i, T_v) for the tilting complex T; for a linear Nakayama module category Hom between indecomposables has dim <= 1 (thin modules), so the question is how the mutation rule can raise dim to 2, 3 (a bound by mutation depth?). That would turn items 3-4 into a theorem. Also: is "b2 grounds at most one of two b1-equal paths" forced by T being a complex over an interval-module category.
- Experimentalist: max dim e_iAe_v by mutation depth at n = 8, 9 (does it grow?), and the dim-3 rows (3 at c0): are they in the 22-type shape? Classes 2-3, n = 9 prefix.
- Skeptic: refute item 3 by searching for any (row, i) with dim >= 3 and a circuit; scalar != 1 pairs (the coefficient-2 relations) are the place a counterexample would live.
- Toolsmith: a hand-built algebra with an nn 2-cycle whose Coxeter key IS an LNA key (the only way D-type can be on a walk): can one exist? If not, that is a proof candidate.
