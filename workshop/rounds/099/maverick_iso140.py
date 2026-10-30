"""Round 099 maverick: DECIDE whether C_a and C_b are Z-congruent (P C_a P^T = C_b, P in GL_n(Z)) for LNAs whose Phi has exactly one Salem pair and S = C + C^T has signature (n-1, 1).
Reduction: C = S (1 - Phi^-1)^-1 (det(1-Phi) != 0), so P C_a P^T = C_b  <=>  R = P^T satisfies R Phi_b = Phi_a R and R^T S_a R = S_b.  R lies in the commutant lattice Lam
(integer R with R Phi_b = Phi_a R, rank 14 here).  The only indefinite part of S is the Salem hyperbolic plane (eigenvectors p, q of Phi, S(p,q) = 1), so the isometries are a finite set times Phi^Z.
Majorant H = S + (S d)(S d)^T, d = q - p.  For a solution R: tr(H_b^-1 R^T H_a R) = tr(H_b^-1 S_b) + z^T S_b H_b^-1 S_b z, z = t q_b - p_b/t (R^-1 p_a = p_b / t),
and multiplying R by a power of Phi_a moves t into [lam^-1/2, lam^1/2]; this gives BOUND.  We enumerate ALL R in Lam with norm <= BOUND + slack (pari qfminim) and test R^T S_a R = S_b exactly.
usage (repo root):  timeout 10m .venv/bin/python workshop/rounds/099/maverick_iso140.py [slack]            (the open pair; also the 140-group sample table)
                    timeout 10m .venv/bin/python workshop/rounds/099/maverick_iso140.py table [per_quipu] [slack]"""
import sys, itertools, collections; sys.path.insert(0, "workshop/rounds/097")
from maverick_cartanZ import *
import numpy as np
N = 12
pari.allocatemem(3 * 10 ** 9)
npf = lambda M: np.array(M.tolist(), dtype=float)
def pq(S, P):
    w, V = np.linalg.eig(P)
    i = [k for k in range(N) if abs(w[k].imag) < 1e-9 and w[k].real > 1.5][0]; j = [k for k in range(N) if abs(w[k].imag) < 1e-9 and abs(w[k].real - 1 / w[i].real) < 1e-9][0]
    p, q = V[:, i].real, V[:, j].real; c = p @ S @ q
    return p, q / c, w[i].real
def major(S, p, q):
    v = S @ (q - p); return S + np.outer(v, v)
def test(a, b, slack=0.5, SC=10 ** 6, verbose=False):
    Ca, Cb = cartan(N, a), cartan(N, b); Pa, Pb = phi(Ca), phi(Cb); Sa, Sb = Ca + Ca.T, Cb + Cb.T
    assert sig(Sa) == (N - 1, 1) and sig(Sb) == (N - 1, 1) and int(Sa.det()) == int(Sb.det())
    Sa_f, Sb_f, Pa_f, Pb_f = map(npf, (Sa, Sb, Pa, Pb))
    pa, qa, lam = pq(Sa_f, Pa_f); pb, qb, lam2 = pq(Sb_f, Pb_f)
    assert abs(lam - lam2) < 1e-9
    Ha, Hb = major(Sa_f, pa, qa), major(Sb_f, pb, qb); Hbi = np.linalg.inv(Hb)
    base = np.trace(Hbi @ Sb_f)
    ts = np.exp(np.linspace(-0.5 * np.log(lam), 0.5 * np.log(lam), 2001))
    BOUND = base + max((t * qb - pb / t) @ Sb_f @ Hbi @ Sb_f @ (t * qb - pb / t) for t in ts) + slack
    n = N; m = n * n
    A = sp.kronecker_product(sp.eye(n), Pb.T) - sp.kronecker_product(Pa, sp.eye(n))   # row-major vec: vec(R Pb - Pa R)
    K = pari.matkerint(pari.matrix(m, m, [int(v) for v in A])); k = int(pari.matsize(K)[1])
    Bs = [np.array([int(K[i, j]) for i in range(m)], dtype=object).reshape(n, n) for j in range(k)]
    Bf = [B.astype(float) for B in Bs]
    G = np.array([[np.trace(Hbi @ Bi.T @ Ha @ Bj) for Bj in Bf] for Bi in Bf])
    Gi = pari.matrix(k, k, [int(round(G[i, j] * SC)) for i in range(k) for j in range(k)])
    U = pari.qflllgram(Gi); Gr = U.mattranspose() * Gi * U
    mind = min(int(Gr[i, i]) for i in range(k)) / SC
    res = pari.qfminim(Gr, int(BOUND * SC) + 10 ** 5, 300000, 0)
    nv = int(pari.matsize(res[2])[1]); assert 2 * nv == int(res[0]) and nv < 150000
    Ui = np.array([[int(U[i, j]) for j in range(k)] for i in range(k)], dtype=object)
    Sa_i = np.array(Sa.tolist(), dtype=object); Sb_i = np.array(Sb.tolist(), dtype=object); hits = []
    for c in range(nv):
        x = Ui.dot(np.array([int(res[2][i, c]) for i in range(k)], dtype=object))
        R = sum((int(x[j]) * Bs[j] for j in range(k)), np.zeros((n, n), dtype=object))
        if (R.T.dot(Sa_i).dot(R) == Sb_i).all(): hits.append(R)
    ok = False
    for R in hits:
        P = R.T; Cai = np.array(Ca.tolist(), dtype=object); Cbi = np.array(Cb.tolist(), dtype=object)
        d = int(pari.matdet(pari.matrix(n, n, [int(v) for v in P.flatten()])))
        if (P.dot(Cai).dot(P.T) == Cbi).all() and abs(d) == 1: ok = True; Pgood = P
    return dict(bound=BOUND, mindiag=mind, rank=k, nvec=nv, nsol=len(hits), congruent=ok, P=(Pgood.tolist() if ok else None))
if __name__ == "__main__":
    if sys.argv[1:2] == ["table"]:
        per = int(sys.argv[2]) if len(sys.argv) > 2 else 3; slack = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
        keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
        for rl, kk in keys.items(): groups[kk].append(rl)
        a0 = (2, 2, 2, 3, 0, 2, 3, 0, 3, 0)
        grp = [v for v in groups.values() if tuple(a0) in map(tuple, v)][0]; print("group size", len(grp))
        q = collections.defaultdict(list)
        for rl in grp:
            r = quipuForAlmostSeparateLNA(N, list(rl))
            if r is not None: q[r].append(tuple(rl))
        print("quipus", {k_: len(v) for k_, v in q.items()})
        reps = [(r, rl) for r in q for rl in q[r][:per]]
        for (r1, a), (r2, b) in itertools.combinations(reps, 2):
            res = test(a, b, slack); print("within" if r1 == r2 else "ACROSS", a, b, {k_: v for k_, v in res.items() if k_ != "P"}, flush=True)
    else:
        slack = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
        a = (2, 2, 2, 3, 0, 2, 3, 0, 3, 0); b = (2, 2, 3, 0, 3, 0, 2, 2, 3, 0)
        print("open pair", test(a, b, slack)); print("control a,a", {k_: v for k_, v in test(a, a, slack).items() if k_ != "P"})
