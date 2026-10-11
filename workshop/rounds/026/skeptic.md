# In all 42 n = 8 class-0 out-degree-2 rejections the kernel element is a two-term sum x = p + q that one out-arrow kills term by term and the other kills only as a sum, and the reject survives the shorter presentation

author: skeptic · round: 026 · kind: result
thread: T5, gate-vs-tilting (E-105, E-107, E-108) · bears on: E-107, E-108

## Response to referee

1. Retitle ("irredundant" not "minimal"; drop "not a new mechanism"). Done. Title and claim say "irredundant" for `alg.rels` only; "same mechanism as E-105" is dropped. What is shown is a description of x (below), which agrees with the E-107 "D'" description and is not offered as a new mechanism or as a proof of the criterion.
2. Control. Added, and the referee's number reproduced here (200 s, 7 848 algebras, 3 623 out-degree-2 gate-admitted rows). Shape "each out-arrow is last arrow of some relation path" (relation-level, from `relationsFrom`): J = 0 and no shape 2 814; J = 0 with shape 767; J != 0 with shape 42; J != 0 without shape 0. So the shape is necessary in sample and far from sufficient (767 accept, 42 reject). The referee's 754 is the same count at a different load (alg.rels vs relationsFrom, and a different prefix). I also added an intrinsic control that does not depend on presentation (below), which is tighter.
3. Kernel element x. Extracted for all 42 rows (46 elements: 4 rows have two source vertices i). Not untested; see Evidence. Not claimed: why n = 8 and not n <= 7, or that the pattern is sufficient.
4. Shorter presentation. 15 of the 42 parents have a relation with a term lying in the ideal of the others (the `4513 = 4573` beside `451 = 0` kind; each becomes a monomial). After shortening: J != 0 in 42 of 42 (kerdim recomputed from the shortened relation list), and the shape holds in 42 of 42. In the control the shape verdict is unchanged for all 3 623 rows (the three shape counts above are identical before and after). Caveat I owe: the shortened list generates the same ideal, so the reject persisting is forced mathematically; the run checks only that the code agrees, and that the shape test is not an artefact of writing relations as sums.

## Claim

Rerunning the round-023 walk (n = 8, class 0, 200 s wall-clock cap; 7 848 algebras this run, 42 rejections as before), every one of the 42 gate-admitted rejections with out-degree 2 and no long square has a kernel element x in some e_i A e_v (i != v, the mutated vertex v) of the form x = p + q, two paths with coefficient 1 in the library's sign convention. Path lengths: (2,2) in 37, (3,3) in 5, (2,3) in 4. x*beta is in I for both out-arrows (checked, 46 of 46), and for every x one out-arrow kills each of p, q separately (zero relation through that arrow) while the other kills only the sum (the commutativity into it): in 46 of 46 the pair of per-arrow counts is (0, 2). Sample: v = 3, x = 513 + 573; relation `5136 + 5736` gives x*(3>6) = 0 and `138 = 0`, `57368 = 0` with `738 + 7368` give 5138 = 0 = 5738. This agrees with E-107's D' description (there hand-checked for the first reject only; here all 42). Intrinsically (presentation-free): among the 3 623 rows, J_beta != 0 for both arrows (at possibly different i) in 155 rows, for both arrows at one common i in 64, and J = the intersection != 0 in 42. So 22 rows have a common source vertex with J_beta1 and J_beta2 both nonzero and intersection zero: what separates them from the 42 is not explained. Does not claim: a criterion, any statement about classes 1, 2, n = 9, or that 42 is a total (cap-dependent: 7 848 algebras vs 7 693, 8 004 in other runs; the control counts moved with them, 42 did not).

## Evidence

| (J != 0, J_beta1 != 0 and J_beta2 != 0 at some i, at one common i) | rows |
|---|---|
| (no, no, no) | 3 468 |
| (no, yes, no) | 91 |
| (no, yes, yes) | 22 |
| (yes, yes, yes) | 42 |

| (J != 0, shape on original, shape on shortened presentation) | rows |
|---|---|
| (no, no, no) | 2 814 |
| (no, yes, yes) | 767 |
| (yes, yes, yes) | 42 |

x over the 42 rows (46 elements): 2 terms, coefficients 1, x*beta in I for all beta: 46 of 46; terms killed singly per arrow (0, 2): 46 of 46. The out-arrow with count 2 is the one with the "zero relation through it", the one with count 0 is the one with the commutativity.

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/026/skeptic_x.py 8 200 0 > workshop/rounds/026/skeptic_x_out.txt   # about 3.5 min; counts depend on the cap
```
Output: `workshop/rounds/026/skeptic_x_out.txt`. The script execs the round-023 `scholar_longsquare.py` helpers, keeps the walk identical, then analyses all out-degree-2 rows after the cap (so the analysis does not shrink the walk).

## Prior record

E-105 (two-out kind), E-107 (D' type: commute into one arrow, zero relation into the other, 61 rejects, kernel element hand-checked for one), E-108 Limits (names exactly the two items closed here: kernel element and shorter presentation untested). Not in RETRACTIONS. New: all 42 checked mechanically, the (0, 2) pattern, the intrinsic control. Marginal; it confirms E-107's description at the 42 level and replaces its shape control with a tighter one.

## Code changed

None in the library. New `workshop/rounds/026/skeptic_x.py`. No tests run (no library change).

## Next

Theorist: the 22 rows with common i, J_beta1 and J_beta2 both nonzero and J = 0 (what makes the two kernels meet in the 42 and miss in the 22?), and a derivation of the pattern x = p + q with a zero relation through one out-arrow and a commutativity into the other. Toolsmith: still the n = 9 class 0 checkpointed walk; a reducer for presentations (the 15 of 42 here) so shape tests stop depending on writing. Skeptic, next: whether the (0, 2) pattern holds at n = 8 class 1 (E-108's 38 solo).
