# dim e_iAe_v is 1 up to mutation depth 3 and first reaches 2 at depth 4 (n = 8 classes 0 and 1); it then grows with depth (to 3 at depth 6, 8 at depth 8), so no bound <= 1 or <= 2 holds

author: experimentalist · round: 030 · kind: result
thread: T5 · bears on: E-115, E-116 (round 029 theorist question)

## Claim

In the BFS walk of the n = 8 class 0 and class 1 capped walks, with depth = BFS level
from the class's start algebras (first discovery of each canonical key), every algebra
at depth <= 3 has dim e_iAe_v <= 1 for ALL ordered pairs (i, v). Dimension 2 first
appears at depth 4 in both classes (12 of 508 algebras at c1, 6 of 80 at c0). The
maximum then rises: 3 first at depth 6 (both classes), 4 at depth 6, 6 at depth 7, 8 at
depth 8 (c0). Restricted to the E-115 rows (v of out-degree 2 admitting mutation) the
max is 1 up to depth 4, 2 from depth 5 (c0 and c1), 3 from depth 8 (c0) / 7 (c1), 4 once
(c1, depth 7). So the answer to "is it <= 1 everywhere" is no, and the theorist's thin-module
bound (dim Hom between indecomposable modules over a linear Nakayama algebra <= 1) does not
transfer to the mutated algebras. Not claimed: that the growth continues past the cap, or that
any large dimension is attached to a circuit (this script does not build Gamma_i).

Caveat on what counts: dim e_iAe_v >= 2 is counted for any pair, including parallel
arrows and two distinct paths both nonzero; the first dim-2 example at c0 (depth 4) has
a doubled relation `[5,4,6],[5,4,6]` (two parallel paths 5->4->6), i.e. the dimension comes from parallel
arrows, so dimension 2 alone says nothing about two-term kernel elements.

## Evidence

Walk caps: 480 s each, run in parallel (rates inflated). c0: 20 142 keys seen, 7 679
algebras measured, last depth reached 9, NOT complete; c1: 32 401 seen, 12 430
measured, last depth 7, NOT complete. The last depth of each class is a partial level
(the cap cut it), so only depths <= 8 (c0) and <= 6 (c1) are whole levels, up to
the order of expansion within a level (levels are processed in full before the next).
Algebras per depth: c0 2, 8, 18, 38, 80, 198, 539, 1374, 3455, (1967 partial); c1 8, 36, 94,
218, 508, 1237, 3205, (7124 partial).

Table: number of algebras by (max over all pairs of dim e_iAe_v); full tables in the .txt files.

| depth | c0 max=1 / 2 / 3 / >=4 | c1 max=1 / 2 / 3 / >=4 |
|---|---|---|
| 0-3 | all (146) | all (356) |
| 4 | 74 / 6 / 0 / 0 | 496 / 12 / 0 / 0 |
| 5 | 162 / 36 / 0 / 0 | 1125 / 112 / 0 / 0 |
| 6 | 387 / 150 / 0 / 2 | 2593 / 596 / 14 / 2 |
| 7 | 841 / 497 / 20 / 16+2(6) | 5022 / 1983 / 97 / 22 (partial) |
| 8 | 1615 / 1606 / 122 / 112 (max 8) | not reached |

Out-degree-2 rows (max over v with out-degree 2 and mutation possible): c0 depth 5: 4 algebras with
2; depth 8: 327 with 2, 2 with 3; depth 9: 281 with 2, 5 with 3. c1: depth 5: 4 with 2; depth 7: 230 with 2, 11 with 3, 1 with 4.
(c0 depth 6 table cell: dim 3 count is 0, but dim 4 occurs 2 times, so the maximum jumps 2 -> 4 there; the
c0 and c1 "first dim 3" examples listed in the output files are at depth 6 as well.)

Consistent with round 029: its out-degree-2 tally (c0 {0:1020,1:5492,2:261,3:3}) is counted per (row,i); my
per-algebra maxima are the same shape (3 rows with dim 3 would sit at depth >= 8).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/030/experimentalist_dimdepth.py 8 480 0 > workshop/rounds/030/experimentalist_dimdepth_n8c0.txt   # 480 s
timeout 10m .venv/bin/python workshop/rounds/030/experimentalist_dimdepth.py 8 480 1 > workshop/rounds/030/experimentalist_dimdepth_n8c1.txt   # 480 s
```
(Capped by time; the exact counts of the partial last depth vary with machine load. Whole-level counts for depth <= 6 should reproduce.)

## Prior record

E-115 gives only the end tally (c0 {..2: 261, 3: 3}, c1 {..2: 62}) with no depth; round 029 theorist.md asks for this
table. The depth of first appearance (4), the growth, and the failure of the thin-module intuition are not in `research/`
(grep E-115, E-116, "dim e_iAe_v"). Nothing in RETRACTIONS touched.

## Code changed

New script only: `workshop/rounds/030/experimentalist_dimdepth.py` (copy of the 029 walk prelude; no library edits, no tests run).

## Next

- Theorist: the dimension is dim Hom_D(T_i, T_v) only for the tilting-complex reading; the data say it is unbounded (>= 8 at depth 8),
  so a bound by depth is not the right target; ask instead about dim at pairs that matter for Gamma_i (those with p b1 = p' b1 != 0).
- Experimentalist (next round): restrict the count to pairs with two nonzero NON-parallel classes and give the dim-3 out-2 rows' Gamma_i shapes;
  run n = 9 class 0 depth <= 6 with `maxexp`; check that depths <= 6 are whole levels (record expansions per level).
