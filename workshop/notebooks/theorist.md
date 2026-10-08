# Theorist notebook (rewritten round 042)

## What I believe now
- Setting: C_B = C' + E_{vi} (J = e_i, E-136; checked on every gate-admitted record at n = 6 (74 384) and n = 7 (44 761)). When column v of C' is e_v (H1) and the v-row of C' is u = e_w - e_i (H2), then C' = [[Z,0],[u^T,1]] with Z = C_A on V\v, and Q(x) = x (adj S_ii - adj S_wi - adj S_iw), S = xZ + Z^T. Exact (P1). H1+H2 hold on 716/766 (n = 6 walk), 82/82 (n = 6 guard off), 120/120 (n = 7).
- With F = Z Z^-T and c_k = chi(S_w, F^k S_w) (c_0 = 1, c_-k = c_k-1), Q_1 and Q_2 are the first two orbit moments of S_w. If F e_w = e_i then Q_1 = 0 for free (triangularity) and Q_2 = 1 + c_2; so Q_2 = 1 iff c_2 = 0, which holds when F e_i = -e_m (proved) and in all 296 s = 1 samples (observed). Other samples have F^-2 e_w = e_i (s = -2, 540 at n = 6): moments c_1 = 0, c_2 = 1, c_3 = 0 observed, unproved.
- Not a matrix identity of Z: unitriangular Z with Y_wi = 1 give B_1 anywhere in -37..13. The law needs realisability (F-orbit of a simple being short, Phi S_w = S_i or Phi^2 S_i = S_w) and the one-arrow shape of the new vertex.
- The n = 4 counterexample (E-141) has u_i = 0, u = e_2 + e_4, outside H2; Q = 0 there. R = 1 corresponds to B = 0.
- Earlier (still standing): J_i = Hom(S_v, e_iA) (E-122); gate = no single nonzero path in J (E-114); dim J_i <= d_i - 1 (E-126); L2 conditional on AI 2.31; x = 0, infinity and trace routes give no obstruction (041).

## What I tried
- 042: `rounds/042/theorist_{dump_off,reduce,lemma,moments,orbits,exceptions,cbcheck,n4,zrandom,explore,blocks}.py`. 041: `rounds/041/theorist_{dump,analyse,diffpoly,d1,matrix,random,example,bfs4}.py`. 039: d3walk, guard, out2. 037: d3. 034: dimji. 031: t3.

## Next
- Prove H2 and the orbit relation from the algebra: build eAe at n = 6 for the s = 1 shapes, find the single relation i ~> v -> w, check Ext^2(S_i, S_w) and Phi S_w = S_i by AR theory (module level, not Cartan level).
- Explain s = -2 (540 of 716 at n = 6): is it the dual of s = 1 with w, i swapped? Test with Z^T and the x^{n-2} end.
- Search for a step with H1, H2 and c_2 != 0 (n = 8, or random realisable Z): that would be the key-preserving candidate (needs Q_2 = 0).
- The 50 off-shape steps (u = e_a + e_b - e_i, or H1 fails) also have Q = x^2 + ...: generalise P1 with m = C_B-row.
- Blind spots: 766 is a time-limited prefix of one walk (314 distinct (Z,w,i)); "dim J = 1 always" is a sample property; types counted from c_k for |s| <= 6 only.
