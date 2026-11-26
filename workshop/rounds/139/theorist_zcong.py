"""Round 139 theorist: independent rebuild of the F-055 / E-269 Cartan-level test (Phi-equivariant Z-congruence).
No import of the round 097/099/110 scripts and no pari.  Own: Cartan builder, Phi, commutant lattice (sympy nullspace +
Smith saturation), majorant, LLL (Gram form) and Fincke-Pohst enumeration, exact integer check; own q_S histogram.

Question: is there R in GL_n(Z) with R Phi_b = Phi_a R and R^T C_a R = C_b  (P = R^T gives P C_a P^T = C_b)?
Bound (re-derived): S = C + C^T has signature (n-1,1); p, q eigenvectors of Phi for lam, 1/lam (isotropic), S(p,q) = 1;
W = <p,q>^perp is S-positive definite.  H = S + (S(q-p))(S(q-p))^T is positive definite, {p, q, S-ON basis of W} is H-ON.
An intertwining isometry sends p_b -> s p_a, q_b -> q_a/s, W_b -> W_a isometrically, so
Q(R) = tr(H_b^-1 R^T H_a R) = s^2 + s^-2 + n - 2;  replacing R by Phi_a^k R (still a solution) puts |s| in
[lam^-1/2, lam^1/2], so some solution has Q(R) <= n - 2 + lam + 1/lam.  Enumerate all of Lam with Q <= bound + slack.

usage (repo root):  timeout 10m .venv/bin/python -W ignore workshop/rounds/139/theorist_zcong.py [slack] [N ...]
"""
import sys, time, random, collections, math
from fractions import Fraction
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_decomp

# the 19 cross pairs of E-292 (copied from workshop/rounds/110/skeptic_qsno.out and skeptic_qsno12.out; group-of-140 row excluded)
PAIRS = {
 12: [((2,2,3,0,3,0,3,0,0,0),(2,4,0,0,4,0,0,2,3,0)), ((2,2,2,2,3,0,3,0,3,0),(2,3,0,2,5,0,0,0,3,0)),
      ((2,3,0,2,3,0,4,0,0,0),(2,4,0,0,2,4,0,0,3,0)), ((2,2,3,0,3,0,2,3,0,0),(2,3,0,4,0,0,2,2,3,0)),
      ((2,2,3,0,2,3,0,3,0,0),(2,2,4,0,0,2,3,0,3,0)), ((2,2,2,2,3,0,3,0,0,0),(2,2,5,0,0,0,2,3,0,0)),
      ((2,2,2,2,2,3,0,3,0,0),(2,2,2,3,0,2,4,0,0,0)), ((2,2,2,2,2,4,0,0,3,0),(2,2,3,0,2,4,0,0,0,0)),
      ((2,2,2,3,0,2,3,0,0,0),(2,2,5,0,0,0,2,2,3,0)), ((2,2,2,2,3,0,2,3,0,0),(2,2,2,4,0,0,2,2,3,0)),
      ((2,2,3,0,2,2,3,0,0,0),(2,2,4,0,0,2,2,2,3,0)), ((2,2,3,0,2,2,2,3,0,0),(2,4,0,0,2,2,2,2,3,0))],
 11: [((2,2,2,3,0,3,0,3,0),(2,3,0,4,0,0,2,3,0)), ((2,2,3,0,2,3,0,3,0),(3,0,2,2,4,0,0,3,0)),
      ((2,2,2,2,2,3,0,3,0),(2,2,4,0,0,2,2,3,0)), ((2,2,2,2,3,0,3,0,0),(2,2,4,0,0,2,3,0,0))],
 10: [((2,2,2,4,0,0,3,0),(2,2,3,0,3,0,0,0)), ((2,2,2,2,3,0,3,0),(2,2,3,0,2,3,0,0))],
 9:  [((2,2,2,3,0,3,0),(2,4,0,0,2,3,0))],
}

def cartan(n, rl):
    """C[j][i] = 1 iff the path i -> j (i <= j) of the linear A_n quiver is nonzero; a relation of length l at vertex s
    (1-based, rl[s-1] = l) kills every path containing the segment s -> s+l."""
    rel = [(s + 1, l) for s, l in enumerate(rl) if l]
    C = [[0] * n for _ in range(n)]
    for i in range(1, n + 1):
        for j in range(i, n + 1):
            if not any(i <= s and s + l <= j for s, l in rel): C[j - 1][i - 1] = 1
    return sp.Matrix(C)

def phi(C): return -(C.T.inv()) * C

def signature(S):
    w = np.linalg.eigvalsh(np.array(S.tolist(), dtype=float)); return (int((w > 1e-9).sum()), int((w < -1e-9).sum()))

def qS_hist(S):
    """Discriminant form of the even lattice S: invariant factors and histogram of q(x) = x^T S x / 2 mod 1 over S^-1 Z^n / Z^n."""
    n = S.shape[0]; D, U, V = smith_normal_decomp(S, domain=sp.ZZ)   # D = U S V
    d = [abs(int(D[i, i])) for i in range(n)]; assert all(x != 0 for x in d)
    nz = [i for i in range(n) if d[i] > 1]
    M = V.T * S * V   # x = V D^-1 z  =>  x^T S x = z^T D^-1 M D^-1 z
    Mf = [[Fraction(int(M[i, j]), d[i] * d[j]) for j in nz] for i in nz]
    h = collections.Counter()
    def rec(k, z):
        if k == len(nz):
            v = sum(z[i] * Mf[i][j] * z[j] for i in range(len(nz)) for j in range(len(nz))); h[(v / 2) % 1] += 1; return
        for t in range(d[nz[k]]): rec(k + 1, z + [t])
    rec(0, [])
    return tuple(sorted(x for x in d if x > 1)), tuple(sorted(h.items()))

def commutant(Pa, Pb):
    """Z-basis of {R integer n x n : R Pb = Pa R} (saturated), as list of integer sympy matrices."""
    n = Pa.shape[0]
    rows = []
    for i in range(n):
        for j in range(n):   # entry (i,j) of R Pb - Pa R ; R[k,l] -> index k*n+l
            r = [0] * (n * n)
            for k in range(n): r[i * n + k] += int(Pb[k, j]); r[k * n + j] -= int(Pa[i, k])
            rows.append(r)
    ns = sp.Matrix(rows).nullspace()
    B = sp.Matrix([[x for x in (v * sp.ilcm(*[e.q for e in v])).T] for v in ns])   # k x n^2 integer rows
    k = B.shape[0]; D, U, V = smith_normal_decomp(B, domain=sp.ZZ)   # D = U B V ; row_i(V^-1) = row_i(U B)/d_i
    UB = U * B; basis = []
    for i in range(k):
        di = int(D[i, i]); r = UB.row(i); assert all(int(x) % di == 0 for x in r)
        basis.append([int(x) // di for x in r])
    from sympy.polys.matrices import DomainMatrix   # exact integer LLL (Euclidean) to tame the Smith-form entries
    dm = DomainMatrix([[sp.ZZ(x) for x in r] for r in basis], (k, n * n), sp.ZZ)
    for delta in (sp.QQ(3, 4), sp.QQ(1, 2)):   # sympy's ddm_lll has hit an internal assert on large entries; retry, else keep basis
        try: dm = dm.lll(delta=delta); break
        except AssertionError: pass
    red = dm.to_Matrix()
    return [sp.Matrix(n, n, list(red.row(i))) for i in range(k)]

def eig_pq(S, P):
    w, Vv = np.linalg.eig(P)
    real = [i for i in range(len(w)) if abs(w[i].imag) < 1e-9]
    i = max(real, key=lambda t: abs(w[t].real)); lam = w[i].real
    j = [t for t in real if abs(w[t].real - 1 / lam) < 1e-9]; assert len(j) == 1 and lam > 1 + 1e-9
    assert sum(1 for t in range(len(w)) if abs(abs(w[t]) - abs(lam)) < 1e-7) == 1   # lam simple
    p, q = Vv[:, i].real, Vv[:, j[0]].real; c = p @ S @ q
    return p, q / c, lam

def lll_gram(G, delta=0.99):
    """LLL on a basis given by its float Gram matrix; returns (unimodular integer T, reduced Gram T^T G T)."""
    k = G.shape[0]; T = np.eye(k, dtype=object); G = G.copy()
    def gso(G):
        mu = np.zeros((k, k)); Bn = np.zeros(k)
        for i in range(k):
            for j in range(i):
                mu[i, j] = (G[i, j] - sum(mu[j, l] * mu[i, l] * Bn[l] for l in range(j))) / Bn[j]
            Bn[i] = G[i, i] - sum(mu[i, l] ** 2 * Bn[l] for l in range(i))
        return mu, Bn
    i = 1
    while i < k:
        mu, Bn = gso(G)
        for j in range(i - 1, -1, -1):
            r = int(round(mu[i, j]))
            if r:
                E = np.eye(k, dtype=object); E[j, i] = -r          # b_i <- b_i - r b_j
                Ef = E.astype(float); G = Ef.T @ G @ Ef; T = T.dot(E); mu, Bn = gso(G)
        if Bn[i] >= (delta - mu[i, i - 1] ** 2) * Bn[i - 1]: i += 1
        else:
            E = np.eye(k, dtype=object); E[[i, i - 1]] = E[[i - 1, i]]
            Ef = E.astype(float); G = Ef.T @ G @ Ef; T = T.dot(E); i = max(i - 1, 1)
    return T, G

def short_vectors(G, bound, cap=10 ** 6):
    """All integer x != 0 (one of +-x) with x^T G x <= bound (Fincke-Pohst on Cholesky)."""
    k = G.shape[0]; L = np.linalg.cholesky(G); Q = np.zeros((k, k))   # G = L L^T ; x^T G x = sum_i qii (x_i + sum_{j>i} qij x_j)^2
    for i in range(k):
        Q[i, i] = L[i, i] ** 2
        for j in range(i + 1, k): Q[i, j] = L[j, i] / L[i, i]
    out = []; x = [0] * k
    def rec(i, rem):
        c = -sum(Q[i, j] * x[j] for j in range(i + 1, k)); r = math.sqrt(max(rem, 0) / Q[i, i])
        for t in range(math.ceil(c - r - 1e-9), math.floor(c + r + 1e-9) + 1):
            x[i] = t; used = Q[i, i] * (t - c) ** 2
            if used > rem + 1e-9: continue
            if i == 0:
                if any(x): out.append(list(x))
            else: rec(i - 1, rem - used)
            if len(out) > cap: raise RuntimeError("cap")
        x[i] = 0
    rec(k - 1, bound)
    seen, res = set(), []
    for v in out:
        t = tuple(v); m = tuple(-a for a in v)
        if m not in seen: seen.add(t); res.append(v)
    return res

def zcong(n, Ca, Cb, slack=0.5):
    Pa, Pb = phi(Ca), phi(Cb); Sa, Sb = Ca + Ca.T, Cb + Cb.T
    assert Ca.det() == 1 and Cb.det() == 1
    assert signature(Sa) == (n - 1, 1) == signature(Sb), "signature"
    if Pa.charpoly() != Pb.charpoly(): return dict(cong=False, why="charpoly", nvec=0, nsol=0, rank=0)
    f = lambda M: np.array(M.tolist(), dtype=float)
    pa, qa, lam = eig_pq(f(Sa), f(Pa)); pb, qb, lam2 = eig_pq(f(Sb), f(Pb)); assert abs(lam - lam2) < 1e-9
    Ha = f(Sa) + np.outer(f(Sa) @ (qa - pa), f(Sa) @ (qa - pa)); Hb = f(Sb) + np.outer(f(Sb) @ (qb - pb), f(Sb) @ (qb - pb))
    Hbi = np.linalg.inv(Hb)
    basis = commutant(Pa, Pb); k = len(basis); Bf = [f(B) for B in basis]
    G = np.array([[np.trace(Hbi @ Bi.T @ Ha @ Bj) for Bj in Bf] for Bi in Bf]); G = (G + G.T) / 2
    bound = n - 2 + lam + 1 / lam + slack
    T, _ = lll_gram(G)
    Bi = [np.array(B.tolist(), dtype=object) for B in basis]
    Bt = [sum((int(T[i, j]) * Bi[i] for i in range(k)), np.zeros((n, n), dtype=object)) for j in range(k)]   # exact new basis
    Btf = [B.astype(float) for B in Bt]
    Gr = np.array([[np.trace(Hbi @ X.T @ Ha @ Y) for Y in Btf] for X in Btf]); Gr = (Gr + Gr.T) / 2
    vecs = short_vectors(Gr, bound)
    Cai, Cbi = np.array(Ca.tolist(), dtype=object), np.array(Cb.tolist(), dtype=object)
    sols = []
    for v in vecs:
        R = sum((int(v[j]) * Bt[j] for j in range(k)), np.zeros((n, n), dtype=object))
        assert (R.dot(np.array(Pb.tolist(), dtype=object)) == np.array(Pa.tolist(), dtype=object).dot(R)).all()
        if (R.T.dot(Cai).dot(R) == Cbi).all(): sols.append(R)
    return dict(cong=bool(sols), nvec=len(vecs), nsol=len(sols), rank=k, lam=lam, bound=bound, R=sols[0] if sols else None)

def controls(n, a):
    """E-292's control: first two members of a's quipu in lnaKeys order (library classification, not the test)."""
    from quivermutation import coxeterTables as ct
    from quivermutation.quipuForms import quipuForAlmostSeparateLNA
    keys = ct.lnaKeys(n); key = keys[tuple(a)]
    qa = quipuForAlmostSeparateLNA(n, list(a)); mem = []
    for rl, kk in keys.items():
        if kk == key and quipuForAlmostSeparateLNA(n, list(rl)) == qa: mem.append(tuple(rl))
        if len(mem) == 2: break
    assert mem[0] == tuple(a), (mem, a)
    return mem[1]

if __name__ == "__main__":
    from quivermutation import coxeterTables as ct
    slack = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    Ns = [int(x) for x in sys.argv[2:]] or [9, 10, 11, 12]
    random.seed(139)
    for n in Ns:
        for a, b in PAIRS[n]:
            t0 = time.time()
            Ca, Cb = cartan(n, a), cartan(n, b)
            lib = all(sp.Matrix(ct.lnaCartanMatrix(n, x)) == cartan(n, x) for x in (a, b))
            ha, hb = qS_hist(Ca + Ca.T), qS_hist(Cb + Cb.T)
            r = zcong(n, Ca, Cb, slack)
            a2 = controls(n, a); C2 = cartan(n, a2); rc = zcong(n, Ca, C2, slack)
            # planted YES: b' = U C_a U^T with U a random product of elementary matrices (true answer YES)
            rp, fails = None, 0   # a U whose conjugate makes the float Gram not positive definite counts as a failed try
            while rp is None and fails < 10:
                U = sp.eye(n)
                for _ in range(4):
                    i, j = random.sample(range(n), 2); E = sp.eye(n); E[i, j] = random.choice([-1, 1]); U = U * E
                try: rp = zcong(n, Ca, U * Ca * U.T, slack)
                except np.linalg.LinAlgError: fails += 1
            rp = rp or dict(cong="undecided", nsol=-1)
            print(f"n={n} a={''.join(map(str,a))} b={''.join(map(str,b))} lam={r['lam']:.4f} |detS|={abs(int((Ca+Ca.T).det()))} "
                  f"libC={lib} qS_differ={ha != hb} cross: cong={r['cong']} nvec={r['nvec']} nsol={r['nsol']} rank={r['rank']} | "
                  f"ctrl a2={''.join(map(str,a2))} cong={rc['cong']} nvec={rc['nvec']} nsol={rc['nsol']} | planted cong={rp['cong']} nsol={rp['nsol']} failedU={fails} "
                  f"{time.time()-t0:.0f}s", flush=True)
