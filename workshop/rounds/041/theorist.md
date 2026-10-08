# "J_i != 0 forces R(x) != 1 + t" is false as a general statement (4-vertex counterexample); it holds only on the LNA-class walks, where the Coxeter polynomials differ by x^2 + ... with leading coefficient 1

author: theorist · round: 041 · kind: negative
thread: T5 · bears on: H-015, E-136, E-138, E-137, E-132

## Claim

Precise statements, with what is proved and what is only computed.

**(P, proved).** Let A be a finite-dimensional path algebra of an acyclic quiver modulo relations, v a gate-admitted vertex whose mutation is legal, B the reduced child with the same vertex set and acyclic. Put C' = r C_A r^T, J_i = dim Hom(S_v, e_iA), H = e_v J^T, so C_B = C' + H (E-136, taken as given). Then det C_A = det C_B = det C' = 1, so t = j^T C'^{-1} e_v = 0 and the key guard passes iff R(x) := det(xC_B + C_B^T)/det(xC' + C'^T) is identically 1. For j = e_i (J supported at one vertex with dim 1), with N = (xC' + C'^T)^{-1}: R(x) - 1 = x N_{iv} + N_{vi} - x (N_{ii} N_{vv} - N_{iv} N_{vi}). Consequences, all by direct expansion: R(0) = R(infinity) = 1 + t = 1, and R(x) = R(1/x). So **x = 0 and x = infinity give no obstruction**, and the notebook's suggested "x^{n-1} coefficient" obstruction is not one: the x^1 (equivalently x^{n-1}) coefficient of det(xC_B + C_B^T) - det(xC' + C'^T) equals det(C')(G_{vi} - (G C'^T G)_{iv} - G_ii G_vv), G = C'^{-1} (formula checked numerically on 156 steps), and it is **0 on all 156 walk steps with J = e_i** (n = 6: 45, n = 7: 111).

**(N, refuted, hand-checkable).** The universal statement "gate-admitted and J != 0 imply the child's Coxeter polynomial differs from the parent's" is **false**. Counterexample, n = 4: arrows 1->2, 1->3, 2->3, 2->4, 3->4, one relation (1-2-3-4) = (1-3-4), v = 3. Gate admits; J_1 = 1 (d_1 = 2; both paths 1->3 composed with 3->4 are identified); child B: arrows 1->2, 2->4 (two arrows), 4->3 with relation (second 2->4 arrow)(4->3) = 0. Both have det C = 1 and Coxeter polynomial x^4 - x^3 - 3x^2 - x + 1 (`theorist_example.py`; hand check: C_A and C_B printed, both rows of paths listed above). So R(x) = 1 there with H != 0. The true statement is therefore not about the matrix identity and cannot be a theorem of the form "R != 1 + t for every gate-admitted J != 0 step".

**(M, matrix level, computed).** No proof from C' alone can exist: for random unitriangular C (entries 0..3, n = 4..6), r = the mutation matrix with S = minimal support of row v, and the realistic constraints d_i = C[i,v] >= 2 and rank bound sum_{w in S} C[i,w] >= d_i - 1, there are C with C' + E_{vi} having the same polynomial as C' (n = 4: 6 of 3 192 trials; n = 5: 67 of 99 914; n = 6: 6 of 57 454; n = 7: 0 of 26 984). These matrices are not claimed to be Cartan matrices of algebras; they only show that the argument must use realisability.

**(W, computed on walks).** On the guarded class-0 walks, every step with J != 0 has J = e_i of dimension 1 at n = 6 (45 steps) and n = 7 (111 steps; 4145 and 10198 expansions), and Q(x) := det(xC_B + C_B^T) - det(xC' + C'^T) is x^2 (1 + 2x + x^2) or x^2 (1 + x + x^2) at n = 6 (23 / 22 steps) and x^2 (1 + x + x^2 + x^3) at n = 7 (all 111): lowest degree exactly 2, coefficient 1. This is a stronger and more informative form of E-138's "the key fails" (the polynomials differ by a very small, regular amount), and it is a walk observation: not proved, n = 8 not run, and n = 6 is the only length with two shapes.

**(L, computed, bears on the key guard).** Random acyclic algebras (not walk parents): of the gate-admitted legal J != 0 steps, those with the child's key equal to the parent's key all have a parent key that is **not an LNA key** (seed 3: n = 4: 9 of 213 J != 0 steps; n = 5: 24 of 516; n = 6: 19 of 522, none on an LNA key; the 59 n = 6 steps from LNA-key parents all change the key; n = 5 seed 1: 5 of 188, parents not LNA-keyed). So the key guard's refusal of J != 0 steps is exact in this sample for parents whose key is an LNA key, but not for arbitrary parents. Not claimed: any statement for n >= 7 off walks, or that this survives non-random parents.

What is **not** claimed: that the walk law ("J != 0 => key moves") is a theorem; the evidence for it remains E-138 (n = 6, 7, 8 class 0) plus (L); no proof was found.

## Evidence

| quantity | n = 6 | n = 7 |
|---|---|---|
| guarded class-0 walk records (gate-admitted legal steps) | 15 273 | 39 836 |
| C_B = r C_A r^T + H holds | 15 273 / 15 273 | not run |
| steps with J != 0 (all J = e_i) | 45 | 111 |
| lowest degree of Q | 2 (45) | 2 (111) |
| dP_1 = 0 and formula agrees | 45 / 45 | 111 / 111 |
| hypothetical C' + E_{vi} or C' + E_{iv} keeping the polynomial, at all 12 971 distinct (C', v) of n = 6 | 0 of 129 710 | not run |

The last row is the strongest "realisable" evidence of the walk law at matrix level: on the walks' own (C', v), changing any single entry (v,i) or (i,v) never keeps the polynomial, although (M) shows that arbitrary C' would. So the obstruction is a property of the Cartan matrices that occur on LNA-class walks, not of the shape of the perturbation.

Random off-walk search (`theorist_random.py`, 150 s per n, seed 3; n = 4 also seed 2, n = 5 also seed 1): tallies in the Reproduction outputs. Child equality uses `search._coxeterKeyOrNone` on parent and child (computed on the child). The weakest step: random relations may be non-minimal presentations (e.g. the n = 5 hit with (1-4-5) = (1-2-4-5)); the n = 4 counterexample is minimal and was checked by hand against the printed Cartan matrices.

## Reproduction

```
.venv/bin/python workshop/rounds/041/theorist_example.py                       # n = 4 counterexample, 2 s
timeout 10m .venv/bin/python workshop/rounds/041/theorist_dump.py 6 0 100 workshop/rounds/041/theorist_dump_n6c0.pkl   # 100 s
timeout 10m .venv/bin/python workshop/rounds/041/theorist_dump.py 7 0 450 workshop/rounds/041/theorist_dump_n7c0.pkl   # 451 s
.venv/bin/python workshop/rounds/041/theorist_analyse.py workshop/rounds/041/theorist_dump_n6c0.pkl   # C_B = C'+H, hypotheticals, 40 s
.venv/bin/python workshop/rounds/041/theorist_diffpoly.py workshop/rounds/041/theorist_dump_n6c0.pkl  # Q(x), 1 s (also n7c0.pkl)
.venv/bin/python workshop/rounds/041/theorist_d1.py workshop/rounds/041/theorist_dump_n6c0.pkl        # x^1 coefficient formula
CON=2 .venv/bin/python workshop/rounds/041/theorist_matrix.py 5 3 60000        # matrix-level, tens of s
timeout 10m .venv/bin/python workshop/rounds/041/theorist_random.py 6 150 3   # 150 s per n
```

## Prior record

E-136 gives C_B = C' + H; E-138 gives the walk law (n = 6, 7, 8 c0) and the R(x) = 1 + t restatement; E-132/E-134/E-137 already record gate-admitted J != 0 steps with the key kept, but at off-walk fans (n = 6) and without a minimal example. New here: the n = 4 minimal example, the (M) negative (no matrix-level proof exists), the Q(x) regularity on walks, the dP_1 = 0 correction, and (L). No retraction touches these (grepped R-numbers via "coxeter"/"preserv").

## Code changed

None in the library. New scripts in `workshop/rounds/041/`: `theorist_dump.py`, `theorist_analyse.py`, `theorist_diffpoly.py`, `theorist_d1.py`, `theorist_matrix.py`, `theorist_random.py`, `theorist_example.py`, `theorist_bfs4.py` (tilting-only BFS from the n = 4 example; did not finish, no result, so whether A and B are tilting-equivalent is open). No tests run (no library file touched).

## Next

- Theorist (or skeptic): explain Q_2 = 1 on walks: it is the next-order obstruction; a proof would show the x^2 coefficient of Q is the invariant "number of ... " (guess: a count of simple-socle elements; unproved). This is the real target; the R(x)-identity route via x = 0, infinity, trace is closed.
- Skeptic: check the n = 4 counterexample independently (build B by hand; Cartan matrices), and whether the 9 + 24 + 18 same-key random hits are all non-LNA keys for a structural reason (my guess: their polynomials are not products of cyclotomics with the LNA shape; untested).
- Experimentalist: if the key-guard-off tables at n = 6..8 contain a J != 0 child with an LNA key, compare its Q(x) to (W); a Q with lowest degree > 2 or coefficient 0 is the candidate for a passing step.
