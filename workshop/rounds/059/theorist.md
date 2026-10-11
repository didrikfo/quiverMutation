# dim End(T) = Cartan(rewrite) at every loopless step checked (8 152, J = 0 or not); End(T) = H_0 and J lives in H_1, so End(T) cannot see J (End(T) = rewrite's algebra conditional on G2)

author: theorist · round: 059 · kind: proof
thread: T10 / T5 · bears on: H-015, E-168, E-174, E-177, E-179
scope: proof part: any finite-dimensional A = kQ/I and loopless v, conditional on Oppermann Thm 1.1 (statement read in the arXiv LaTeX, proof not checked) and on one identification (G2 below). Experiment: all LNAs and duals n = 4..8, every vertex with an out-arrow, no gate (8 152 steps (30 + 112 + 420 + 1584 + 6006), no parallel arrows, acyclic); 38 steps with 7 parallel-arrow parents (n = 7, random unguarded walk depth <= 8, 600-node frontier cap, seed 59). Dimension level only (Cartan matrix of End(T) vs Cartan(child), uniformly transposed); the algebra-level statement is conditional on G2. Not covered: cyclic quivers, loops, n >= 9, class-2 or n = 10 edges.

## Response to referee

1. Stale output: regenerated `theorist_enddim_out.txt` from the current script (n=6, 7, 8; 0 mismatches each; n=8 2m35). Labels now read `End=child(H^T only)` on every step, as the text says. The script also prints a new SPLIT line.
2. Claim: "Hence" is now conditional on G2 (Lemma item (3) and the sentence after it); title and experiment sentences say "dimension". Proved/conditional: End(T) = H_0 (Oppermann Thm 1.1, statement only); empirical: G2, supported by dim only.
3. Sign convention: in the script `res[m]` = dim Hom(T,T[m]). Hom(T,T[1]) = 0 at all 8 152 steps (SPLIT line: `Hp==0` everywhere; this is silting, AI convention). J_i = Hom(T,T[-1]) (`Hm`), nonzero exactly on the J != 0 steps (2574 at n=8). `J0` in the script requires Hm = 0 and Hp = 0; Hp is redundant, so J0 <=> Hm = 0. Item (1) below is stated in that convention.
4. `D = -Hm` test: it accepts H or H^T independently of eqH/eqHt. New SPLIT counts (n=6,7,8, all steps): on the J = 0 steps (252, 924, 3432) D=-Hm holds for both H and H^T (Hm = 0, so no information); on J != 0 steps (168, 660, 2574) it holds for H^T only, never for H, and eqHt holds at all steps while eqH holds at none. So the uniform transposition is supported by both tests exactly on the J != 0 steps, and by eqH/eqHt alone elsewhere.
5. Optional (hand-checked A5 G2 instance): not done this sitting; needs a hand computation of Oppermann degree-2 generators, left for next round.
Also: Thm 1.1 quote is from the LaTeX as served, not stored in `research/literature/sources/` (I may not edit `research/`); re-fetch command is under Reproduction.

## Claim

**Lemma.** Let A = kQ/I be finite dimensional, v a vertex with no loop, T = T_v + sum_{j != v} P_j with T_v = (P_v -> sum_{alpha: v->h} P_h) in degrees (-1, 0) (left AI mutation). Let (kM, d_M) be the dg quiver of Oppermann Section 3 built from the minimal dg model of A. Then
(1) T is silting (E-168; AI 2.31). Convention: Hom(T,T[m]) is computed as in round 050; Hom(T,T[1]) = 0 (silting; observed 0 at every step) and Hom(T,T[m]) = 0 for m <= -2 (T lives in two adjacent degrees).
(2) REnd(T) ~ (kM, d_M), and (kM, d_M) is homologically non-negative. So End_K(T) = H^0 = H_0(kM, d_M) = kM_0 / (d_M of the degree-1 arrows), and Hom(T,T[-1]) = H_1.
(3) (CONDITIONAL, = G2, unproved) The repo's steps 1-7 plus `reduce` present H_0(kM, d_M).
Verified: **dim End(T) = Cartan(rewrite) (uniformly transposed) at all 8 152 + 38 steps, J = 0 or not.** Given (3), **End(T) = the rewrite's algebra for every loopless v, with or without J = 0**; without (3) only (1), (2) and the dimension check stand. J_i = Hom(S_v, e_i A) (AI 2.50(c)) = Hom(T,T[-1])_{(v,i)} enters only through H_1: T is tilting iff H_1 = 0 (all J_i = 0), and then REnd(T) is formal, equal to End(T), derived equivalent to A. J != 0 changes H_1 and nothing in H_0.
**Corollary (explains E-174).** The End(T) comparison (dimension or iso) cannot discriminate J = 0 from J != 0, as a theorem, not as an accident of 25 steps. J = 0 must be tested by Hom(T,T[-1]) or J_i directly.
For acyclic A: Cartan(End T) = R C R^T + Hom(T,T[-1]) entrywise, because the Euler form is R C R^T, the negative Homs of m <= -2 vanish and the positive ones vanish; Hom(T,T[-1]) sits in the (v, i) entries (dim J_i), as E-138 has it (see Prior record).
It does **not** claim: that the rewrite equals H_0 on cyclic quivers (step 4 of Oppermann is not implemented, F-002); that the child is derived equivalent to the parent when J != 0 (it is not); that anything about *chains* of steps follows beyond "each step is End of the one-step T of its own parent".

## Evidence

**Proof sketch, weakest steps marked.**
(1) E-168. Hom(T,T[<=-2]) = 0 is degree reasons: T in K^[-1,0], T[-m] in K^[m-1,m].
(2) Theorem 1.1 of arXiv:1504.02617 (LaTeX fetched this round, `mutation.tex`, version as served 2015-04-10, ~v1/v2 UNVERIFIED which; statement matched): "(kQ,d) dg quiver algebra, homologically non-negative, i a vertex with no loops of degree 0. Then the derived endomorphism ring of the left Aihara-Iyama silting mutation of kQ at i is (kM, d)." Hypotheses: only no degree-0 loop at i; parallel arrows, relations of any shape, higher generators are all allowed. Non-negativity of M: arrows into i go up, arrows out of i go down but have degree >= 1 unless in A (these become alpha* of degree 0), composites have degree ||phi|| - 1 = |phi| >= 0, anti-composites have degree >= 1. H_0 of a non-negative dg quiver with d lowering degree is kM_0 / (d of kM_1), where the degree-1 elements are generated by the degree-1 arrows: nothing else enters.
(3) **G2, weakest step.** Dictionary (E-lit 1504.02617 summary): step 1 = composites alpha phi, step 2 = flips alpha*, step 3 = r-bar (a relation out of i is a degree-1 arrow that rotates to degree 0), step 4 = relation sum alpha* alpha beta, step 5 = r-bar alpha* = r/alpha is d(phi alpha^-1), step 6 = relations into i. Step 7 (the kernel over all paths out of i*) is claimed to be what remains: the d-images of the degree-2 generators x (relations among relations = overlaps, Oppermann Ex. 4.2) that start at i and become degree 1 in M. **I did not prove that step 7 equals those d(x).** Only the data below bear on it.
(4) Conventions: left versus right (the repo's `mutateLeftAtVertex` dualises) and Oppermann's composition order (alpha phi = phi first). The experiment accepts H or H^T; re-run at n = 6 splitting them: Cartan(child) = H^T at every step and H at none that is not also symmetric-equal, one uniform transposition (the repo's Cartan[i][j] versus the round-050 Hom(T_i,T_j)), never a mix. A uniform transpose is what a left/right (op) convention difference would give; it is not a mismatch of algebras but is not excluded as an op.

**Experiment A (`theorist_enddim.py`).** For every LNA and dual of length N and every vertex v with an out-arrow, no gate: child = `procedure.mutateAtVertex` (steps 1-7, not reduced), Cartan(child) = `ap.cartanMatrix`; independently H_ij = dim Hom_K(T_i, T_j[0]) from the round-050 chain-map code (rank mod 2^31 - 1), with R C R^T from `cartanDiscrepancy`. Result: Cartan(child) = H (or H^T) at every step, and R C R^T - Cartan(child) = -dim Hom(T,T[-1]) at every step.

| N | steps J = 0 | steps J != 0 | End = child | mismatches |
|---|---|---|---|---|
| 4 | 20 | 10 | 30 | 0 |
| 5 | 70 | 42 | 112 | 0 |
| 6 | 252 | 168 | 420 | 0 |
| 7 | 924 | 660 | 1584 | 0 |
| 8 | 3432 | 2574 | 6006 | 0 |

Here J = 0 means Hom(T,T[1]) = Hom(T,T[-1]) = 0 in the same code. On J != 0 steps the Euler congruence fails (E-087, E-095) while Cartan(child) = dim End(T) still holds: this is the lemma in its sharp form. The count of relation-set completeness: the child's relation set is built by steps 1-7, so a missing relation would enlarge Cartan(child) over dim End(T); none did (answers E-177's "relation set and dim are inputs" at the dimension level, for these parents). Parity of lna/dual counts is a duality symmetry of the enumeration, not tested for independence.

**Experiment B (`theorist_enddim2.py 7 8 600`, ~4 min).** Parents that are reduced children of unguarded walks, with a parallel pair; 7 parents, 38 steps (28 J = 0, 10 J != 0, 8 with out-degree > 1): End = child, D = -Hom(T,T[-1]), 0 mismatches. Thin: E-174 has 56 parallel-arrow edges plus 8 parallel-arrow failing steps at the iso level; this adds only a dimension-level check on other parents.

**What a parallel-arrow step needs.** Not a new hypothesis in the lemma. Needs: (a) A = the set of arrows out of v *individually*, with composites alpha beta one per pair (the repo's arrow keys do this since 2026-09-18); (b) the dg model uses a *basis* of minimal relations per (source, target) pair, so a relation that mixes parallel paths is one generator, and End(T) fixes the algebra only up to the GL_m action on each block of m parallel arrows, so "equals" means "isomorphic via matrix-valued identification" (as E-174 did), never equality of coefficients; (c) `reduce` must substitute a degree-0 term of d(omega) as an arrow, which for parallel composites needs the coefficient solved, not guessed; (d) no loop at v after earlier steps: parallel arrows do not create loops, but oriented cycles would (step 4, not implemented).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/059/theorist_enddim.py 7     # 30 s
timeout 10m .venv/bin/python workshop/rounds/059/theorist_enddim.py 8     # 167 s
timeout 10m .venv/bin/python workshop/rounds/059/theorist_enddim2.py 7 8 600   # ~4 min
```
Output of the first two: `workshop/rounds/059/theorist_enddim_out.txt`; the third: `theorist_enddim2_out.txt`. Source for Thm 1.1: `curl -sL https://arxiv.org/e-print/1504.02617` (gzip of `mutation.tex`; not committed to `research/literature/sources/`).

## Prior record

E-168 (generation, per step), E-174 (End(T) iso at 25 failing steps + 370 + 35 edges; "Hom(T,T[-1]) != 0 so End cannot see J": observed), E-177 (dim not asserted), E-179 (T5: End(T) = rewrite on J = 0, open), E-138 and E-095 (the Cartan defect is a single (v,i) entry = dim J_i, observed on guarded walks, "no independent End(T)"). **New here:** (i) the reason End(T) cannot see J (H_0 vs H_1) given Oppermann Thm 1.1, with its hypotheses checked against the LaTeX; the lemma is therefore stronger than the J = 0 statement E-179 asked for, and the J = 0 case needs nothing extra; (ii) an independent End(T) dimension check on every vertex of every LNA n <= 8 without the guard (extends E-138's 1-entry defect from guarded walks to all vertices and includes the independent Hom computation E-138 lacked). Not new: the Cartan defect formula itself (E-138).

## Code changed

None. New scripts `workshop/rounds/059/theorist_enddim.py`, `theorist_enddim2.py` (import round-050 `skeptic_tilt.py` code by exec). No tests touched.

## Next

- scholar: check Oppermann Section 3 and Example 4.2 against `procedure` step 7 (G2): does the kernel in step 7 equal d of the degree-2 generators out of i? A hand example: A5 with two overlapping relations, mutate at the start of the overlap.
- toolsmith: `compare2` dim assertion (E-177) is now supported by this dimension identity; a version that also runs on parallel-arrow edges of the n = 7 c1 and n = 10 E-169 paths would close dim-level completeness there.
- skeptic: the match is with H^T at all steps (one uniform op); decide whether that is the intended left/right convention of `tiltingPlus` or a bug in the reading of Cartan.
- Open: cyclic quivers (Oppermann step 4) are outside everything above; chains of steps need nothing new if each parent is an algebra, but a parent of a J != 0 step is not derived equivalent to the start.
