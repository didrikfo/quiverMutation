# Every J_i != 0 row of E-137 is a dead-end step (the child leaves the Coxeter class), and out(i) = 3 is not forced by the gate: it is the smallest shape of a "born" J_i with d >= 3, not a theorem

author: theorist · round: 039 · kind: negative (with one exact criterion)
thread: T5 · bears on: E-131, E-133, E-134, E-137, E-138, H-015

## Claim

(1) **Re-reading the table.** E-131/E-133/E-137 tabulate (d_i, J_i) at every gate-admitted v of every algebra the key-preserving walk *visits*, whether or not the child of the step at v stays in the class. Re-running E-137's walk with that one extra column (n = 8, c0, 7 852 expansions, 574 s): all 229 rows with J_i != 0 (192 steps (A, v)), including all 9 rows with d >= 3 (rows 7798, 7810, 7822, 7831, 7836 reproduce), have a child whose key differs from the class key. All 121 rows with d >= 3 and J_i = 0 pass. So the rows with J_i != 0 are not on the walk; they are the edges the walk refuses. E-131's "d_i = 2 on walks" and E-137's "out(i) = 3" describe parents of refused steps, not a property of the class. Same at n = 6 and 7 (class 0, 100 s and 200 s): 65 and 58 steps with H != 0, 0 pass.

(2) **out(i) = 3 has no gate-level reason.** Hand algebra (theorist_out2.py): vertices 1..7, arrows 1->2, 1->3, 2->4, 2->5, 4->6, 5->6, 3->6, 6->7; relation 1-2-4-6-7 = 1-3-6-7 (the path 1-2-5-6-7 free). It is gate-admitted at v = 6 with d_1 = 3, dim J_1 = 1, out(1) = 2 (and with a second commutation dim J_1 = 2). Its key is not an LNA key and its child key differs from its own. So "out(i) = 3 whenever J_i != 0, d_i >= 3" is false for the gate and is a statement about this sample.

(3) **Why the sample has out(i) = 3 (mechanism, heuristic, not proved).** Let i have out-arrows a_1..a_m. Then e_iAe_v = sum_k a_k e_{t(a_k)}Ae_v, and an element x = sum a_k y_k of J_i = {x : x * rad = 0} is either inherited (if only one a_k has y_k b != 0, then y_k lies in J_{t(a_k)}, or a relation through a_k kills it) or born at i, which needs a relation between at least two branches (a commutation, E-133's p1 b = p2 b). In every d >= 3 row (9 of 9) J_i is born at i (no out-neighbour has J != 0 at the same v); the d = 2 rows are either inherited (out 1: 24 rows) or squares (out 2: 176 born). A born element with d >= 3 needs three paths in e_iAe_v: either three branches of one path each (out 3, T1 of E-134) or two branches one of which has two paths (my hand example, out 2, 7 vertices). The first is smaller, and the walk at n = 8 reaches only algebras with dim A 47..75 in its last two BFS levels. Row 7798 by hand (not machine-verified): v = 7, i = 5, out(5) = {1, 3, 4}, paths 5-1-7, 5-3-7, 5-4-1-7, 5-4-3-7 with 5-1-7 = -5-3-7, so d = 3, and the kernel element is p(517) + p(5417) + p(5437): it uses all three out-arrows of 5. Not claimed: that out(i) = 3 fails at some walk; only that no argument forces it.

(4) **When does a J != 0 child pass the key guard (E-138).** Write C' = r C_A r^T (congruent to C_A, so chi_{C'} = chi_A), j = (dim J_i)_i, e = e_v, C_B = C' + e j^T, M(x) = x C' + C'^T, N = M^{-1}, chi_X(x) = det(x C_X + C_X^T)/det C_X. Matrix determinant lemma (rank 2: x e j^T + j e^T = [e j][x j e]^T) gives, with a = j^T N e, f(x) = x a(x), b = j^T N j, c = e^T N e, and e-term = e^T N j = f(1/x) (because M(x)^T = x M(1/x)):

  det M_B / det M' = (1 + f(x)) (1 + f(1/x)) - x b(x) c(x) =: R(x),  det C_B / det C' = 1 + t, t = j^T C'^{-1} e.

**The key guard passes iff R(x) = 1 + t identically.** Necessary: t = 0 for derived equivalence (det C_B = det C_A). Exact check on all 123 H != 0 steps at n = 6, 7: the lemma reproduces det(xC_B + C_B^T)/det(xC' + C'^T) at x = 2, 3, 5 (tolerance 1e-6), and det C_B = det C_A = 1 in every one (t = 0), so the determinant alone never rules a step out; the polynomial does. Also t = -(r j)^T C_A^{-1} e (r is an involution), so t = 0 reads sum_w j_w (C_A^{-1})_{wv} = 0 (convention of C_A^{-1} not rechecked). Not claimed: that no J != 0 child ever passes. 0 of 352 observed steps pass; I have no theorem and no example either way.

## Evidence

Tables (data, not counted). n = 8 c0, rows saved when J_i != 0 or d >= 3:

| (d, J != 0?) | child passes key guard | rows |
|---|---|---|
| d >= 3, J = 0 | yes | 121 |
| d >= 3, J != 0 | no | 9 |
| d = 2, J != 0 | no | 220 |

(d, out(i), J_i born or inherited at the same (A, v)), J != 0 rows: (2,1,born) 9, (2,1,inh) 24, (2,2,born) 176, (2,2,inh) 7, (2,3,born) 4, (3,3,born) 2, (4,3,born) 4, (5,3,born) 3.
Guard at n = 6 / n = 7: (guard passed, det C_A, det C_B, sum dim J, lemma ok): (False, 1, 1, 1, True) with counts 65 / 58.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/039/theorist_d3walk.py 8 0 570 12000 workshop/rounds/039/theorist_d3walk_n8c0.pkl   # 574 s, load-dependent
.venv/bin/python workshop/rounds/039/theorist_d3rows.py workshop/rounds/039/theorist_d3walk_n8c0.pkl   # 2 s; output in theorist_d3rows.txt
.venv/bin/python workshop/rounds/039/theorist_out2.py                                                  # 10 s
timeout 10m .venv/bin/python workshop/rounds/039/theorist_guard.py 6 0 100   # 100 s;  and  ... 7 0 200  # 200 s
```

## Prior record

E-137 (the out(i) = 3 table: unchanged, I add the guard column), E-138 (formula for C_B; I use it), E-116 (gate admits the two-term kernel), E-131 (d = 2 "on walks"). Grep of EXPERIMENTS.md/FINDINGS.md for the pass/fail of J != 0 children: no record that all J != 0 steps leave the class (E-138 reports 60 H != 0 steps, only the identity). The determinant-lemma form R(x) is not in `research/`. Not in RETRACTIONS. Weakest part: (3) is a plausibility argument; (1) is for one class and one capped walk (the rows are the last two BFS levels), so "never passes" is an observation.

## Code changed

None in the library. New: `theorist_d3walk.py` (E-137's script plus edges, relations, guard flag), `theorist_d3rows.py`, `theorist_out2.py`, `theorist_guard.py` (E-138 referee script plus pass flag, det, lemma check). No tests touched.

## Next

- skeptic: look for a gate-admitted J != 0 step whose child passes the key guard (any n <= 7, parents not necessarily on a walk); one such case would break (1) as a law. Also check the hand-computed kernel element of row 7798 with `skeptic_kernel.py` on this pkl's algebra (edges/rels are saved).
- theorist (me): prove or refute "J != 0 implies R(x) != 1 + t" in the single-i case (j = e_i), where R is a 2 x 2 minor of N; try x = -1 and the x^{n-1} coefficient (a trace condition) as first obstructions.
- experimentalist: rerun E-131's (d, J) table counting only children that pass; the d = 2 vs d >= 3 question then has no J != 0 rows left to explain.
