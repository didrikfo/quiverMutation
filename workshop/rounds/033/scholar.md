# AI 2.32(b) at a vertex is exactly "J_i = Hom(S_v, e_iA) = 0"; derived here, and it matches step 7 of E-066 (J_8 = span of the commutativity element, dim 1); it adds no obstruction to circuits

author: scholar · round: 033 · kind: result (restatement + check, partly rediscovery)
thread: T5 · bears on: H-015, E-066, E-078, E-095, E-121

## Claim

Convention: right modules, P_x = e_xA, paths composed left to right, Hom(P_y,P_x) = e_xAe_y. For a basic finite-dimensional A, a vertex v with no loop, and D = (+)_{j!=v} P_j:
1. The minimal right add(D)-approximation of P_v is g: D' = (+)_{b: v->t(b)} P_{t(b)} -> P_v, left multiplication by the arrows b out of v (the projective cover of rad P_v; any P_j -> P_v, j != v, lands in rad P_v and lifts). Parallel arrows count separately.
2. AI 2.31: N := [D' -g-> P_v] (D' in degree 0, P_v in degree 1; this is cone(g)[-1], the triangle N -> D' -> P_v -> N[1] of Def. 2.30 for mu^-) with D gives a silting object, always.
3. AI 2.32(b) for M = add A: mu^- is tilting iff Hom(g, P_i): Hom(P_v,P_i) -> Hom(D',P_i), x |-> (x b)_b, is injective for every i != v (for M in D the condition is vacuous).
4. ker of that map = {x in e_iAe_v : x b = 0 for all out-arrows b} = {x : x rad A = 0} (other arrows have wrong tail) = Hom(P_v/rad P_v, P_i) = Hom(S_v, e_iA) =: J_i = the S_v-socle of e_iA.
5. Degree sense: J_i = H^{-1} of the complex Hom(N, P_i) = Hom_D(N, P_i[-1]). It is NOT the cohomology H^{-1}(N) (N has H^0 = ker g, H^1 = coker g = S_v). Link: truncation triangle ker g -> N -> S_v[-1]; Hom(S_v[-1], P_i[-1]) = Hom(S_v, P_i) -> Hom(N, P_i[-1]) is an isomorphism because Ext^{-1}, Ext^{-2} of modules vanish. The wording "H^{-1}(cone)" of E-121 should read "H^{-1} RHom(cone, A) = Hom(S_v, A)".
6. Hence: Hom(T,T[m]) for m<0 is nonzero exactly through J (and Hom(D,N[<0]) = 0, since Hom(P_i, g) is onto); silting-not-tilting <=> some J_i != 0. The other piece of the complex, coker_i = Hom(N, P_i) in degree 0, is allowed to be nonzero and is 0 in all computed cases.

Hypotheses: none beyond finite-dimensional, basic, v loopless. Monomial / two-term / scalar-1 is NOT needed for items 1-6 (it is needed only in E-110 to read J != 0 off a circuit graph). So the statement applies to every walk-reachable algebra, not only LNAs. Naming: in this convention the repo's "left mutation at v" (arrows out of v) is AI's mu^- (right approximation, 2.32(b)); the literature file 1009.3370.md also reads 2.32(b) this way. Not claimed: that this explains why walks have no circuits >= 3 (see Prior record).

## Evidence

`scholar_socle.py` recomputes, independently of the `tiltingPlus` call, for every i != v: dim e_iAe_v, kernel of the map over out-arrows (J_map), kernel over all arrows with tail v (J_soc, the socle count), the target dimension, the cokernel, and the Euler identity dim e_iAe_v - dim Hom(D',P_i) = J_i - coker_i. It also prints Cartan(child) - R C R^T.

| parent | v | i | dim e_iAe_v | J_map | J_soc | coker | Euler | Cartan defect |
|---|---|---|---|---|---|---|---|---|
| E-066 step 7 (10 vertices, relation 8,6,4,9 + 8,10,4,9, only out-arrow 4>9) | 4 | 8 | 2 | 1 | 1 | 0 | 1 = 1 | +1 at (4,8) only |
| same, all other i (1,2,3,5,6,7,10) | 4 | - | 0 or 1 | 0 | 0 | 0 | ok | - |
| E-078 n = 5, abde = acde | d | a | 2 | 1 | 1 | 0 | 1 = 1 | +1 at (d,a) only |
| same, i = b, c | d | - | 1 | 0 | 0 | 0 | ok | - |

On step 7 the kernel element is E-066's c = [8,6,4] + [8,10,4]: c * (4>9) = 0 and 4>9 is the only arrow out of 4, so c generates a copy of S_4 in soc(e_8A); that copy is the whole of H^{-1} Hom(N, A). J_map = J_soc and the Euler identity hold at every (parent, i) (assertions in the script). The rewrite's Cartan defect sits exactly in row v at the i with J_i != 0, equal to dim J_i (sign: child minus R C R^T is +1 here; E-095 states X - Y = -dim ker, a different orientation of the difference, not rechecked).

Limits: two parents only (the question's "one explicit mutation" plus a second); the script uses the repo's `idealBasis` for the quotient, so "independent" means independent of `tiltingPlus` and of the rewrite, not of the ideal arithmetic. The derivation of items 1-5 is mine (standard homological algebra, short); it was not compared with the printed AI text (arXiv blocked; I used the repo's literature file, whose item (b) agrees). The least certain step: that the repo's mutation vertex (arrows out of v) corresponds to the approximation of P_v from the right rather than the left in AI's own convention; the computation does not depend on the name.

## Reproduction

```
.venv/bin/python workshop/rounds/033/scholar_socle.py     # about 20 s
```

## Prior record

Recorded: the reduction J_i = {x : x rad A = 0} = Hom(S_v, e_iA) = H^{-1} of the mutation cone, "derived, not tested numerically" (EXPERIMENTS.md E-121 text, line ~12; the AI 2.31/2.32 check "not done"); AI 2.32(b) = Ladkani 2.3(c) = tiltingPlus as one map (E-066, "appears to be, not shown"); kernel location at step 7 (E-066); n = 5 instance (E-078); Cartan defect = -dim ker (E-095). New: the derivation of the one-map identity from 2.32(b) (items 1-4: approximation is the projective cover of rad P_v, so the identification is an equality, not an appearance); the degree correction (item 5); the independence from monomiality; the numerical J_soc = J_map and Euler check. So this closes E-066's "no derivation" limit and the "not done" check, but is mostly a rediscovery as far as consequences go. It does not give a new reason for circuits to be absent: E-121 already concludes the kernel structure cannot forbid nn or circuits >= 3, and the socle reading says only that J_i != 0 is a socle phenomenon of the parent (J intrinsic to the parent, E-097). Whether a socle S_v in e_iA can occur for an algebra derived equivalent to an LNA with a circuit graph of >= 3 edges is exactly the open T5 question.

## Code changed

None to existing files. New: `workshop/rounds/033/scholar_socle.py` (reads `workshop/rounds/001/scholar_h015.py`). No tests run (no repo file touched).

## Next

- theorist: use item 4 to ask, for an LNA-derived parent, which pairs (i, v) can have S_v in soc(e_iA) at all: soc(e_iA) for a derived-LNA algebra is constrained by Hom(S_v, A) = Ext^0(S_v, A); try gl.dim and Hom(S_v, A[m]) = 0 for m != 0 considerations (the tilting condition on the PARENT's own A is automatic; so the obstruction must be on the child's relation to an LNA, still open).
- skeptic: refute item 5 (degree) on a case with coker != 0 (none found among the two parents; a non-gate-admitted vertex with an out-arrow into a vertex with several arrows might give one).
- chair: wording fix for E-121 ("H^{-1}RHom(cone,A)"); naming note on left/right mutation in GLOSSARY.
