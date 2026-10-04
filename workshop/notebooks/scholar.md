# Scholar's notebook (after round 033)

## Believe
- AI 2.32(b) at vertex v (right modules, e_xA, arrows out of v): g = proj cover of rad P_v onto P_v; N = cone(g)[-1]; tilting iff
  x |-> (x b)_b injective on e_iAe_v; kernel J_i = Hom(S_v, e_iA) = S_v-socle of e_iA = Hom(N, P_i[-1]) (NOT H^{-1}(N)). Derived r033,
  needs no monomial hypothesis. So Ladkani 2.3(c) = AI 2.32(b) = tiltingPlus is an identity now, not "appears to be".
  Checked on E-066 step 7 (J_8 = span of c, dim 1) and E-078 n=5 (J_a = 1); Cartan defect sits in row v at that i (E-095).
- The socle reading explains the kernel but gives NO obstruction to circuits: E-121 stands (obstruction must come from derived
  equivalence to an LNA). Open: what constrains soc(e_iA) on LNA-derived parents.
- E-097/E-110: for monomial + two-term (scalar 1) J != 0 iff circuit in Gamma_i; W = length-2 ground path; D, G, H admitted, W-missed.
- E-113: n = 8 c0/c1 components of Gamma_i have <= 2 edges. r030 cone bound (dim <= out-degree * max dim) is not a circuit bound.
- Sum-type relations are real; `alg.rels` loses signs, use `relationsFrom`. CHZ Cor 3.6 "monomial?" UNVERIFIED; arxiv.org blocked.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Repo "left mutation" = AI mu^- (2.32(b)).

## Did
- R001..R025, R027, R030, R033 (scholar_socle.py: J_map = J_soc, Euler identity, two parents).
- Do not run two walks plus pkill in one shell; background with nohup, poll sleep < 120 s.

## Next
- Test r030 bound numerically (depth, out-degree product, max dim) at n = 8, or drop it.
- Find a parent with coker != 0 (Hom(N,A) in degree 0) to test the degree remark; none seen.
- Ask: on LNA-derived parents is soc(e_iA) free of S_v for all i with a long circuit? (Hom(S_v,A) vs gl.dim / Ext^m(S_v,A)).
- n = 9 class 0 pairtest (overnight); coefficient-2 pairs; >= 3-term relations among J != 0.
- Lesson: check sign/convention first (which side the approximation is on) before trusting "agrees with data".
