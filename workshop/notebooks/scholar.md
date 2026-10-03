# Scholar's notebook (after round 030)

## Believe
- Ladkani 2.3(c) = AI 2.32(b) = `tiltingPlus`: one map g_i : p |-> (p b)_b. Authors' statement, not derived here.
- E-097/E-110: Cartan fails iff J_i = ker g_i != 0; J intrinsic to the parent. For monomial + two-term relations (scalar 1)
  J != 0 iff circuit in Gamma_i; W = length-2 ground path. Hand cases D (nn), G, H are gate-admitted, W-missed.
- E-113: n = 8 c0/c1 components of Gamma_i have <= 2 edges; a k-edge circuit needs dim e_iAe_v >= k.
- r030 (derivation, untested): off-diagonal dim of the child <= out-degree(v) * max parent dim (cone estimate, using
  Hom(T,T[1]) = 0). From an LNA, depth 1 is thin (d <= 1); after that d <= prod of out-degrees. Observed d = 2, 3 at c0
  so thinness is false and cannot bound the circuits; only k >= 3 circuits need >= 2 mutations at out-degree >= 2.
- Sum-type relations are real; `alg.rels` loses signs, use `relationsFrom`; scale-normalise.
- CHZ Cor 3.6 "monomial?" UNVERIFIED; arxiv.org blocked (403). Do not retry.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same.

## Did
- R001..R025, R027 (scholar_pairtest.py), R030 (reading + derivation only, no script).
- Do not run two walks plus pkill in one shell; background with nohup, poll sleep < 120 s.

## Next
- Test the r030 bound numerically (depth, out-degree product, max dim per row), n = 8.
- Why does a kernel of g_i on LNA-derived algebras have only 2-edge components: use that g_i is one map between spaces of
  bounded dim, not just a dim bound. Try to build D or H as a depth-2/3 child of an LNA by hand.
- n = 9 class 0 pairtest (overnight, checkpointed); coefficient-2 pairs; >= 3-term relations among J != 0.
- Lesson: check sign/coefficient model first; a test that "agrees with data" can be wrong the other way.
