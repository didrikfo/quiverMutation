# G2 as a span statement: repo step 7 spans red(Z_1 e_i), the degree-2 generators out of i of Oppermann's model (acyclic; 35 + 573 pairs)

author: scholar · round: 061 · kind: result
thread: T5 (T10) · bears on: H-015, E-181, E-179
scope: acyclic loopless quivers, any vertex i. SPAN-level statement. The Oppermann side of the computation is MY model of Z_1 e_i (kernel of d on words p r q in the free path algebra, red = drop q != empty), my reading of Construction 3.1 / Def 3.2, not a recomputation of the paper's output; Ex. 4.2's final answer was NOT recomputed. Computed on 7 hand quivers at n = 5 (35 pairs, all generators single-term) and 573 seeded random pairs at n = 5, 6 (19 with a multi-term generator). Sign/orientation relative to the paper is read, not tested. Not covered: cyclic quivers, parallel arrows, Q^(>=3), and whether Q^(2) of the repo's algebra is a model of the algebra at all (that is Construction 2.2 of the paper, taken as given).

## Response to referee

1. Done by dropping: "same sign"/"no transpose" are out of title and Claim; the script compares spans and cannot see scalars. Sign/orientation is stated as read from the paper (Evidence 1), untested. A scalar-detecting test against Oppermann's explicit d x of Ex. 4.2 at vertex 1 was not done: the paper prints d x only before mutation, so the comparison would again be against my model.
2. Corrected: the table row now says there is no degree-1 cycle out of 1 (1345 is in I; the cause is the missing r_345).
3. Done: seeded random loop folded into `scholar_g2.py` (573 pairs, 19 multi-term generators, 0 differ); the 35 hand pairs are stated single-term. The referee's larger 2 786-pair run is his, not mine.
4. Done in Scope: the Z_1 e_i model is my reading; Ex. 4.2's final answer not recomputed.
5. Not done this sitting (generator count versus minimal count of red(Z_1 e_i) mod left multiples); "generators" in the Claim should be read as "spanning set" until then.

## Claim

Equal as spans. For acyclic Q with admissible I and a vertex i, the relations that step 7 of `procedure._relationsOutOfMutatedVertex` produces (the kernel over the ideal) span, target by target, exactly the set { dec(red(d x)) : x a degree-2 arrow i -> v of Oppermann's model, i.e. a cycle of degree 1 starting at i } , which is the partial of the degree-1 arrows out of i* in M. Conventions, from the paper side only (reading, NOT tested by the script, which compares spans and so cannot see a scalar per generator): `p/alpha` is "paths beginning with alpha" on both sides, the flips alpha*: j -> i* point the same way, and in d(p r q) = (-1)^{|p|} p (d r) q the sign is + for |p| = 0, so by that reading the relative coefficients of `s_P rbar_P` agree with those of `s_P r_P`; the multi-term random cases (19) do test relative coefficients inside a generator against my model. It would be refuted by one acyclic (Q, I, i) where the two spans differ.

Not claimed: that End(T) = the repo's algebra beyond this step (steps 4-6 are matched below only by reading, not by computation), nor anything for J != 0 failures; this closes the "G2" identification in E-181, which remains conditional on Oppermann Thm 1.1 (statement only).

## Evidence

**1. Paper side (read in the LaTeX, `main.tex` of arXiv e-print 1504.02617 = "Quivers for silting mutation"; section/theorem numbers counted by hand from `\newtheorem` (shared counter, `numberwithin` section), not compiled).**
Thm 1.1 (main theorem). Section 3: Construction 3.1 (red, `p/alpha`, dec), Definition 3.2 (the differential). The A = arrows of degree 0 out of i; flips alpha*: j -> i; for phi: i -> j not in A, ||phi|| = |phi| - 1 (+1 if it ends at i); Def 3.2: partial(phi) = dec(red(d phi)) if t(phi) != i. Section 7, proposition "Cancelling anti-compositions": after cancelling phi alpha^{-1} the reduced algebra has arrows phi, alpha phi, alpha* and relations  phi alpha* = -(-1)^{|phi|} (d phi)/alpha. Nothing else (no other arrow kind) enters H_0.
Degree-1 arrows of M, and the repo step that is their partial:
- beta degree 0 into i, ||beta|| = 1: partial beta = sum alpha* (alpha beta)  -> repo step 4 (same sign +);
- phi degree 1 out of i, ||phi|| = 0 (rbar): the relation phi alpha* = +(d phi)/alpha (|phi| = 1) -> step 5 (`rbar alpha* = r/alpha`, same sign);
- alpha phi, phi degree 1 into i, ||alpha phi|| = 1: partial = dec(alpha red(d phi)) -> step 6;
- x degree 2 out of i, ||x|| = 1: partial x = dec(red(d x)) -> **step 7** (this item);
- relations not meeting i: through-relations. (phi alpha^{-1} are cancelled; alpha phi beta* needs a loop at i.)
In Q, d x has degree 1, so it is a sum of words p r q with one degree-1 arrow r (rho_r = d r in kQ_0), p, q degree-0 paths. Words with q nonempty end in some alpha in A (every degree-0 arrow out of i is in A), so `red` keeps exactly the words with q empty: red(d x) = sum_P eps_P s_P r_P, r_P out of i. No Delta occurs (a path out of i that is acyclic never enters i first). Hence partial x = sum_P eps_P s_P rbar_P.

**2. Hand derivation of the kernel.** In the tensor algebra d(p r q) = p rho_r q (|p| = 0). Let z = sum eps_P s_P r_P + sum_alpha y_alpha alpha. Then d z = 0 iff sum_P eps_P s_P rho_P + sum_alpha d(y_alpha) alpha = 0 in the free path algebra. Every path of rho_P starts with some alpha, so this splits alpha by alpha:  sum_P eps_P s_P (rho_P/alpha) = - d(y_alpha). Such y_alpha exists iff the left side lies in the ideal generated by all rho (d of a word a r b is a rho_r b, so image of d on degree-1 elements = I, ker part: I = (rho) by construction). That is the repo's condition "sum eps_P (r_P/alpha) s_P in I for every alpha" (the repo's docstring), word for word. Conversely z in Z_1 with red(z) = 0 is a sum y_alpha alpha with each y_alpha a cycle (alphas are free right factors), i.e. a right multiple of a cycle from another vertex, which the paper's generators of Z_1 supply and which is not out of i*. The repo's `_isForcedByNearer` drops left multiples by paths, the bimodule structure of Z_1. The step I am least sure of: that a minimal set of degree-2 arrows is exactly "Z_1 e_i modulo Z_1 A and modulo left multiples", which uses Construction 2.2 (pick generators of H_1 of Q^(1)); H_1 = Z_1 because Q^(1) has no degree-2 element.

**3. Computation** (`scholar_g2.py`). Oppermann side computed from scratch: kernel of d on all words p r q i -> v in the free path algebra, then red. Repo side: `_relationsOutOfMutatedVertex` fed the same arrow naming `mutateAtVertex` builds, tails read back to old paths (`_oldArrowReading`), extended on the left by every old path to each target v. Compared dimension of K_opp(v), of the repo span, and of their union, per v. Result, 35 hand pairs (every generator single-term) plus 573 seeded random pairs (script loop, seeds 1, 2; 19 with a multi-term generator), all EQUAL (dimensions agree and equal the union):

| Q, I | i with a step-7 generator | where |
|---|---|---|
| A5, 123=0, 234=0 | 1 | v=4 (the generator rbar_{123} followed by 3->4), v=5 forced by it |
| A5, 123, 234, 345 | 1, 2 | v=4,5 ; v=5 |
| A5, 1234, 2345 | 1 | v=5 |
| A5, 123, 2345 | 1 | v=5 |
| Oppermann Ex. 4.2 quiver (1->2->3->4, 2->5), 123, 234 | 1 | v=4 |
| square 1->{2,3}->4->5, (124)-(134), 245, 345 | 1 | v=5 |
| square with (124)-(134), 245 only | none | no degree-1 cycle out of 1: d(12 r_245 - r_c 45) = 1345 (in I, but a boundary only if r_345 existed), and nothing else kills it |

Elements, e.g. A5 with 123, 234 at i=1: repo gives +1 on (r_{123}, 3->4), Oppermann kernel gives -1 on the same basis element, i.e. the same relation `(1* -> 3 -> 4) = 0` up to scalar; in terms of the paper's example, x: 1 -> 4 with d x = s alpha - gamma r (page of Ex. 4.2) goes to partial x = gamma rbar (plus the term ending in alpha, cancelled with s beta^{-1}). **Ex. 4.2 itself mutates at vertex 2, where there is no degree-2 arrow starting at 2: step 7 is empty on both sides** (rows i=2 above), so that example tests steps 4-6 only; the instance that tests step 7 is mutation at vertex 1 of the same quiver. The paper's Ex. 4.2 final answer (partial x = s alpha - s delta^{-1} delta alpha, partial alpha = delta* delta alpha, partial(s delta^{-1}) = s delta*) was not recomputed.

## Reproduction

```
.venv/bin/python workshop/rounds/061/scholar_g2.py      # about 10 s, last line "ALL EQUAL"
```
LaTeX: `curl -sSL https://arxiv.org/e-print/1504.02617 -o e.tgz` (gzip of one file `mutation.tex`; unpacked under the scratchpad, not in `sources/`).

## Prior record

E-181 / E-179 / STATE T5: G2 stated open; hand A5 instance not done. Grep of `research/` for "G2" finds only E-181. Not previously recorded: the identification of step 7 with partial of degree-2 generators, and the reading of steps 4-6 off Definition 3.2 / Section 7 proposition. `research/literature/sources` has no Oppermann LaTeX yet (README row says "seven-step rewrite on a dg quiver"); the chair may promote `main.tex` (78 817 bytes) there. Consequence for T5: the conditional in E-181 loses G2; what remains is Thm 1.1 (statement matched, proof not read), generation by T (E-168), and that End(T) = H_0 cannot see J (E-181), so the live item of T5 is not the identification but the premise J = 0, which this does not touch.

## Code changed

None (new script only: `workshop/rounds/061/scholar_g2.py`; it imports private `procedure._Arrows` and `_relationsOutOfMutatedVertex`).

## Next

- theorist: the one gap in section 2 above (minimal degree-2 generators = Z_1 e_i mod Z_1 A and left multiples) as a short lemma; and whether the Q^(2) truncation is legitimate when H_2(Q^(2)) != 0 (it does not affect H_0 of M, I argued, not proved).
- toolsmith: extend `scholar_g2.py` to the repo's real walk (random n=7 LNA steps, E-151 failing steps) as a regression test; parallel arrows need arrow-keyed words.
- chair: decide whether to put the LaTeX into `research/literature/sources`.
