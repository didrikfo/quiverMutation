# Theorist notebook (rewritten round 041)

## What I believe now
- J_i = Hom(S_v, e_iA) (E-122); gate = "no single nonzero path in J" (E-114); dim J_i <= d_i - 1 (E-126); C_B = r C_A r^T + e_v J^T (E-136). Key guard passes iff R(x) = det(xC_B+C_B^T)/det(xC'+C'^T) is identically 1 (det = 1 throughout, so t = 0).
- Round 041: "gate-admitted J != 0 => R != 1" is FALSE in general. Hand example n = 4: arrows 1->2,1->3,2->3,2->4,3->4, relation 1-2-3-4 = 1-3-4, v = 3: J_1 = 1, child (1->2, 2=>4, 4->3, relation a2.c = 0) has the same Coxeter polynomial x^4-x^3-3x^2-x+1. Random acyclic algebras n = 4..6: ~4% of J != 0 steps keep the key, always on non-LNA keys (0 of ~80 with an LNA-key parent).
- No matrix-level proof exists: R(x) = 1 + x N_iv + N_vi - x(...) has no obstruction at x = 0, infinity, or the x^1 / x^{n-1} (trace) coefficient (that change is 0 on all 156 walk steps). Random unitriangular C' with realistic constraints admit key-preserving perturbations E_{vi}. A proof must use realisability or the walk class.
- Walk law (E-138) refined: on n = 6, 7 class 0, Q(x) = P_B - P' has lowest degree exactly 2 with coefficient 1 (x^2(1+x+x^2+x^3) at n = 7; two shapes at n = 6). Observation only. Under the walk's own (C', v) (12 971 of them at n = 6) no single-entry change keeps the polynomial.
- Earlier: L2 conditional on AI 2.31; 16 key-coincidence fans (E-123) meet an LNA (E-134), unsupported via non-tilting step (E-137).

## What I tried
- 041: `rounds/041/theorist_{dump,analyse,diffpoly,d1,matrix,random,example,bfs4}.py`. 039: `theorist_d3walk.py`, `theorist_guard.py`, `theorist_out2.py`. 037: `theorist_d3.py`. 034: `theorist_dimji.py`. 031: `theorist_t3.py`.

## Next
- Explain Q_2 = 1 on walks (x^2 coefficient of Q): find the algebraic meaning (socle element count? Euler form of the two-term complex?). This is the real target; R(x)-route at x = 0, infinity, trace is closed.
- Off-walk: why do same-key J != 0 hits have non-LNA keys? Test with n = 7 random; check if their polynomials fail the LNA shape (product of cyclotomics?).
- Check `rounds/041/theorist_bfs4.py` result (tilting-only BFS from the n = 4 example A: does B appear?) -- was still running at submission time; if B is reached, A and B are tilting-equivalent and the example is not a "wrong" step in disguise.
- Blind spots: random generator is not uniform and may give non-minimal presentations; the n = 4 example was hand-checked only through printed Cartan matrices; "never on an LNA key" is a sample count.
