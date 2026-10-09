# Cartan-congruence invariants are non-informative about class membership of the 25 key-keeping, tiltingPlus-failing n = 7 children (c1, c2): equal key already implies an integral Cartan congruence P and (here) equal invariants, so the T10 (i) question stays open

author: skeptic · round: 046 · kind: negative
thread: T10 (i) · bears on: H-015, E-149, E-145, STEERING q3
scope: n = 7, key classes 1 and 2 (sizes 12, 14; control: all 66 + 91 LNA pairs congruent, so P existence and invariant equality are non-discriminating at this n; no known out-of-class equal-key pair at n = 7, power of the invariants against the relevant case untested), E-149 walk rebuilt (BFS, key guard, 20 000 expansions each, 459 s / 341 s); all 16 + 9 distinct (parent, v) key-keeping failures of `tiltingPlus` (parents at depth 7-8, all J != 0, |J| = 1 in 23, 2 in 2); 600 J = 0 controls at depth 0-4 only; n = 8, class 0 and the E-145 walk itself not rerun.

## Response to referee

Referee verdict: minor revision. Required items 1-3 done; the claim is weakened accordingly (the test is non-informative about class membership).

1. Control, all LNA pairs (`skeptic_pairs.py 7 1 1 2`, box 1, 3 s). Class 1: 12 distinct Cartan matrices, 66 pairs, P found 66, none 0, timeout 0. Class 2: 14 matrices, 91 pairs, P found 91, none 0, timeout 0. So within an n = 7 key class any two algebras already have an integral P with P C_X P^T = C_Y. P existence for the 25 children is implied by "key kept" and carries no discriminating power; it is no support for membership in the class or for the key guard. The only content of that part is that the step's own map R is not such a P (E-145). Likewise the invariant table (0 of 25 differ) is expected: the LNA key classes at n = 7 carry one signature each (`skeptic_power`), so no separating invariant exists inside them, and equality for the children is what equal key plus this fact predicts. Title and Claim reworded; "every one has an integral P" no longer stands as a result, only as a consequence of equal key.
2. Power. No known out-of-class pair with equal key exists at n = 7 (key-coarser pairs start at n = 12, E-077, and are not proven inequivalent). So the power of the invariants against the relevant pair type is untested at this n, not only that of the back-search. The random-unitriangular power check (`skeptic_randpower`) shows the invariants separate generic forms that are not Cartan matrices of algebras; it says nothing about Cartan matrices with equal key and is not used in the claim. "Do not separate" is a null with no calibrated power.
3. Which c1 steps were checked how. Invariants (`skeptic_inv.py`): all 16 c1 distinct failing key-keepers (0 differ), and all 9 c2 (the referee re-ran c2 only). P search (`skeptic_iso.py`, box 1): c1 15 of 16 found for C_B; c1 step 12 (index 12 in the list) none in box 1 for both C_B and C_B^T, found for C_B in box 2 and 3 (`skeptic_iso12.py`). C_B^T (transposed form): none in box 1 for c1 steps 8, 10 and 12 (a non-find in a box proves nothing; steps 8 and 10 have P for C_B). c2: 9 of 9 found for both C_B and C_B^T. The claim uses C_B only. The 8 c1 parents with parallel arrows were not rebuilt by hand; they rest on the `skeptic_collect` pickle only, as before.

Net: the test cannot say whether the children leave the class. Nothing here is evidence for or against the key guard.

## Claim

At n = 7, classes 1 and 2, none of the 25 distinct key-keeping children that fail `tiltingPlus` is shown to lie outside the derived class of its parent, and this test is non-informative about it. Every congruence invariant tried is equal for child and parent (0 of 25 differ), and for all 25 an explicit P in GL_7(Z) with P C_B P^T = C_A exists (entries in [-1,1] for 24, in [-2,2] for one), although the step's own map R is not one (Cartan incongruent, as in E-145). Both facts are implied by equal key here: all 157 pairs of LNA Cartan matrices in c1 and c2 are congruent (response, point 1). No known out-of-class equal-key pair exists at n = 7, so the invariants' power against the relevant case is untested. The answer to "do they leave the class?" is: not decidable by the tests here; they neither support nor refute it. Since the parents are tilting-only descendants of LNAs (E-149: none tainted), parent and class representative share Cartan congruence class, and so do the children. This does NOT claim the children are derived equivalent to the class: Cartan congruence is necessary, not sufficient, and a bounded tilting-only search from the children failed to reach an LNA, but so did the positive control. A refutation would be one child with an unequal invariant (none found), or a proof that a child is not tilting-reachable.

## Evidence

Invariants (`skeptic_inv.py`), each invariant under every P in GL_n(Z), so inequality proves a different derived class (for finite global dimension), equality proves nothing: Coxeter charpoly; Smith forms of xC + C^T for x = -3..3 (x = 1 symmetrised form, x = -1 skew form); Smith forms of f(Phi), f(Phi)^2 for each irreducible factor f of the charpoly, Phi = -C^-T C; signature of C + C^T; histogram of v^T C v mod m on (Z/m)^7, m = 2..7.

| set | steps | invariant differs |
|---|---|---|
| c1 failing key-keepers (distinct parent, v) | 16 | 0 |
| c2 failing key-keepers | 9 | 0 |
| J = 0 controls (R-congruent by construction), c1 + c2 | 600 | 0 |

Soundness of the invariant set on generic forms only (`skeptic_randpower.py`, random unitriangular 7x7, which are not Cartan matrices of algebras; no claim about Cartan matrices rests on it): 15 of 15 random congruent pairs P C P^T keep all invariants; of 100 pairs with equal charpoly, 23 are separated by another invariant (Smith of xC + C^T most often, q mod 4 next). So the set is not vacuous on generic forms; here it finds nothing. Within LNA key classes at n = 7 (classes 0-5, sizes 8..108, `skeptic_power.py`) all LNAs share one signature, so the key classes carry no further invariant split; E-077 had already found no Smith separation between the key-coarser orbits.

Congruence search (`skeptic_iso.py`, backtracking over vectors f_i with f_i^T C_B f_j = C_A[i,j], entries in a box): P found for 24 of 25 in box 1 (c1 step 12 needs box 2; found also at box 3; C_B^T none in box 1 for c1 steps 8, 10, 12). Control `skeptic_pairs.py`: all LNA pairs, 66 (c1) and 91 (c2), P found in box 1, so P is expected from equal key. Controls: 300 J = 0 steps all found in box 1 (c1 sample, 8 of 8 shown); negative control (random pair, equal charpoly, invariants differing) gives none. Reading: the isometry of the Euler lattices exists, but it is not the map of the step, and nothing here says it is induced by a derived equivalence. Note D = 0 (E-145) is key kept by definition, so this adds no information on D.

Hand rebuild from arrows and relations only (`skeptic_rebuild.py`): c2 9 of 9, c1 8 of 8 rebuildable (the other 8 c1 parents have parallel arrows, which the hand constructor `build` rejects: unchanged limitation from r043); every one: gate True, `tiltingPlus` False, J != 0, key(parent) = key(child), child Cartan matrix equals the recorded one, R C_A R^T != C_B. 300 c2 controls: tiltingPlus True, J = 0, R-congruent. Those 17 failing steps reproduce E-149/E-145 by an independent rebuild; I did not match them one by one with the 13 + 9 distinct E-145 steps (different walk, 500 s guarded versus 20 000 expansions), only the c2 count (9 and 9) and E-149's c1 16 agree.

Bounded tilting-only reach (`skeptic_back.py`): from each rebuilt c2 child (9) and four controls, BFS with J = 0 steps and key kept, 6000 expansions, looking for an LNA or dual of the class: no hit in all 13. The control (children at forward distance 4, so an LNA is at most 4 steps back) also missed, so this test has no demonstrated power and I draw no conclusion from it (consistent with E-142, E-147: reverse search loses about 10% of edges and wanders).

Which invariant can decide: none of these can prove membership. They can prove exclusion, and did not. What could decide membership: an explicit sequence of tilting steps from the child to a class member (long: depth >= 8 both ways), or a derived invariant outside K0 (Hochschild cohomology beyond HH^0, which is trivial for connected acyclic algebras; HH^1 dimension; not in the library).

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/046/skeptic_collect.py 7 1 20000 /tmp/c1.pkl 300   # 459 s; same for class 2 (341 s)
timeout 10m .venv/bin/python workshop/rounds/046/skeptic_inv.py /tmp/c1.pkl                          # 90 s each
timeout 10m .venv/bin/python workshop/rounds/046/skeptic_iso.py /tmp/c1.pkl 1 fail                   # seconds; iso12.py for step 12
timeout 10m .venv/bin/python workshop/rounds/046/skeptic_pairs.py 7 1 1 2                              # control: 66 + 91 pairs, seconds
timeout 10m .venv/bin/python workshop/rounds/046/skeptic_rebuild.py /tmp/c2.pkl                      # fail; add 'ctrl' for controls
timeout 10m .venv/bin/python workshop/rounds/046/skeptic_back.py /tmp/c2.pkl fail 0 4 6000          # 90 s per child
timeout 10m .venv/bin/python workshop/rounds/046/skeptic_randpower.py 400; ... skeptic_power.py 7 6
```
Saved summary: `workshop/rounds/046/skeptic_data.txt` (the pickles are not committed).

## Prior record

E-149 (16 and 9 failures; "not shown: that the failing children leave the derived class"), E-145 (Cartan incongruent under R, key kept, D = 0), E-134/E-137 (non-tilting step is not an equivalence), E-077 (Smith forms do not separate key-coarser orbits), E-142/E-147 (reverse search is weak). The existence of other integral congruences P for all 25 and the equality of the finer invariants are new as statements but follow from equal key (control above); small, non-discriminating.

## Code changed

None in the library. New scripts in `workshop/rounds/046/`: skeptic_collect, _inv, _power, _randpower, _iso, _iso12, _rebuild, _back, _pairs. No tests run (no library file touched).

## Next

For the chair: q3 (promote `tiltingPlus` as a keyword) is not made urgent by these data: the failing key-keepers are not shown to be outside the class, so the key guard is not refuted at n = 7 c1, c2 either; it also is not supported (equal Euler lattices are necessary only). Toolsmith: a reverse tilting search with a positive control at depth >= 7, or an OVERNIGHT job meeting a child with the LNA side at total depth ~16 (E-142 style, with the 10% edge loss fixed). Theorist: is the P found induced by a derived autoequivalence of the child's tilting class (e.g. a composition of other-vertex reflections)? A derived invariant beyond K0 (HH^1 dimension for bound quiver algebras) would be the first test able to separate.
