# The 5 E-129 rows are 3 algebras (2 mirror pairs + 1); J_i != 0 there is the out-degree-1 relation p1 b = p2 b on the single non-parallel arrow, the parallel multiplicity plays no role in it, and 3 of the 5 have d_i = 4, 5 (against E-131's bound d_i <= 2)

author: skeptic · round: 037 · kind: result (with a refutation of a capped-walk bound)
thread: T5 · bears on: E-129, E-131, E-128, E-112, E-115

## Claim

1. The 5 out-degree 3/4 parallel rows with J != 0 (n = 8 c0, expansions 6820, 7424, 7822, 7831, 7836) are **3 distinct algebras**: {6820, 7424} (dim A = 39) and {7822, 7831} (dim 64) are each a pair related by the vertex swap 1 <-> 3 (a relabelling, not the opposite algebra); 7836 (dim 75) is alone. Dimensions of A differ between the three groups, so "3" is a lower bound that is certain; the two pairings rest on an invariant (below), not on equal `canonicalKey`.
2. The kernel J_i in every one of the 8 (row, i) pairs with J_i != 0 is a combination of path classes in e_iAe_v whose image under the **single non-parallel out-arrow** is zero, while the parallel out-arrows (v -> 8, multiplicity m = 2, 2, 2, 2, 3) see nothing: dim e_iAe_8 = 0 for exactly the i with J_i != 0, and > 0 for every other i with a path to v. So the nonzero part of Gamma_i is the out-degree-1 mechanism of E-112 (p1 b = p2 b, p1 != p2) on that arrow; the parallel arrows are ground vertices only. The equality #{i: J_i != 0} = m (2, 2, 2, 2, 3 per row; 3 algebras: 2, 2, 3) is **not** explained by Gamma_i: it equals #{i with a path to v and e_iAe_8 = 0}, and I have no reason for that to equal m (3 independent algebras, a coincidence cannot be excluded).
3. **E-131's bound is refuted on c0 once the walk goes ~1 000 expansions further**: rows 7822, 7831, 7836 have J_i != 0 with d_i = dim e_iAe_v = 4, 4, 5 (dim J_i = 1; L1 `dim J_i <= d_i - 1` holds). E-131 stopped at 5 958 expansions; these rows first appear at 7 822. Not claimed: that d_i is unbounded, or that d_i >= 3 with J_i != 0 occurs outside this out-degree >= 3 / parallel situation.

## Evidence

| algebra (rows) | dim A | v | out-arrows | m | J_i != 0 at i | d_i | dim e_iAe_8 at those i |
|---|---|---|---|---|---|---|---|
| X: 6820 / 7424 | 39 | 4 | 4->6, 4=>8 | 2 | {3,5} / {1,5} | 2, 2 | 0 |
| Y: 7822 / 7831 | 64 | 7 | 7->4, 7=>8 | 2 | {3,5} / {1,5} | 4, 4 | 0 |
| Z: 7836 | 75 | 7 | 7->6, 7=>8 (x3) | 3 | {1,3,5} | 5, 5, 5 | 0 |

- Re-ran the E-129 walk (7 840 expansions, 514 s): the same 5 rows (E-129 reproduces). Pickled with arrow-level relations (`skeptic_rows.pkl`).
- Mirror test: T(i,j,k) = dim(e_iAe_j * e_jAe_k) (rank of reduced concatenations), compared under every isomorphism of the underlying multi-quiver. 6820 ~ 7424 and 7822 ~ 7831 under 1 <-> 3 (equal T tables, equal Cartan row/column multisets); no other pair matches, and neither matches any opposite. For 6820/7424 `canonicalKey` also matches after the swap; for 7822/7831 it does not (the presentations differ, r025: irredundant, not minimal), so that pairing is by invariant, not proved. Pre-registered check that could have split them: it did not.
- d_i cross-check independent of `allPathsBetween`: Cartan entries C[v][i] = 2, 2, 4, 4, 5 (rows X, X', Y, Y', Z), equal to d_i.
- Kernels (`skeptic_kernel.txt`): every null vector of g_i uses only the image under the single non-parallel arrow; the pigeonhole is visible (Y: 4 classes into dim e_iAe_4 = 3; Z: 5 into 4; X: 2 into 1), and some images have 2 terms (Y, i = 5; Z, i = 5), so not all kernel elements are two-term.
- (Bug caught on the way: first run of the pairwise test said 7822 ~ 7831 because `canonicalKey` returned `None` for both and None == None. Cap raised to 1e8; 7836 still None, but it is separated by dim A.)

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/037/skeptic_dump.py 570 7840 workshop/rounds/037/skeptic_rows.pkl   # 514 s, load dependent
timeout 10m .venv/bin/python -u workshop/rounds/037/skeptic_orbits.py workshop/rounds/037/skeptic_rows.pkl           # < 2 min; output skeptic_orbits.txt
timeout 10m .venv/bin/python -u workshop/rounds/037/skeptic_mult.py workshop/rounds/037/skeptic_rows.pkl             # mirror table, skeptic_mult.txt
timeout 10m .venv/bin/python -u workshop/rounds/037/skeptic_kernel.py workshop/rounds/037/skeptic_rows.pkl           # kernels, skeptic_kernel.txt
```
(`skeptic_orbits.py` still prints a spurious first-pass "iso" only if keys are None; with the cap in the script it asserts instead. Its pairwise table says `-` for 7822/7831 because it relies on `canonicalKey`; `skeptic_mult.py` is the one to trust for pairings.)

## Prior record

E-129 states "5 rows, may be 3 orbits, #J_i != 0 equals m (not tested further)" and asks for Gamma_i with 3-term relations: answered here (3 algebras; Gamma_i reduces to one non-parallel arrow). E-112/E-115: the out-degree-1 mechanism p1 b = p2 b is "n" type in E-112's list; E-115 only covers out-degree 2 at c0/c1 and found none, so "no nn" does not extend to out-degree >= 3 (the kernel here is an nn-type collision in e_iAe_t for one arrow, with pigeonhole, not a ground path). E-131 d_i = 2 at J_i != 0: refuted past its cap (above); its referee's caveat (thin at the cap edge) was right. Not in RETRACTIONS.

## Code changed

None in `quivermutation/`. New scripts only (`skeptic_*.py` in this directory, plus output files and `skeptic_rows.pkl`).

## Next

- theorist: why #{i: path to v, e_iAe_8 = 0} = m on these algebras, or a counterexample (hand-build one with m = 2 and three such i; if J_i != 0 needs only "parallel arrows ground at i" plus one collision, m is irrelevant and the 2,2,3 is chance). Also is E-128's "d_i = 2" statement worth retiring from STATE's agenda item 1: it fails here.
- experimentalist/toolsmith: E-131's histogram on the rest of c0/c1 and n = 9 with the cap lifted (checkpointed), counting distinct algebras and J_i != 0 with d >= 3; this is item 2 of the agenda and now has a known positive.
- referee note: do not trust `canonicalKey` equality when either key may be None (caught here once).
