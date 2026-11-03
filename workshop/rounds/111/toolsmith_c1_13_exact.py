"""Round 111 toolsmith, T10: EXHAUSTIVE decision of Z-congruence P C_A P^T = C_B (P in GL_7(Z)) for a key-kept J != 0 step of n = 7 key class 1,
after workshop/rounds/099/maverick_iso140.py (E-269) but with no pari: own integer kernel (saturated commutant lattice of Phi), own LLL and Fincke-Pohst.
S = C + C^T must have signature (6,1) (here: det -4, Phi char poly (x+1)(x^6-x^4-2x^3-x^2+1), one Salem pair).  R = P^T satisfies R Phi_B = Phi_A R and
R^T S_A R = S_B; every solution after multiplying by a power of Phi_A has tr(H_B^-1 R^T H_A R) <= n - 2 + lam + 1/lam (+slack).  Enumerate all R in the commutant
with that norm and test exactly.  Float majorant, not interval arithmetic (as E-269).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/111/toolsmith_c1_13_exact.py [slack] [only_missing]"""
import sys, pickle, collections, random, itertools, time
from fractions import Fraction
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
import numpy as np, sympy as sp
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
N = 7
def phi(C): M = sp.Matrix(C.tolist()); return np.array((-(M.T).inv() * M).tolist(), dtype=object).astype(np.int64)
def int_kernel(A):
    """saturated integer kernel basis (columns x with A x = 0), A integer m x k, by unimodular row reduction of [A^T | I]"""
    m, k = A.shape; M = [[int(A[i, j]) for i in range(m)] + [1 if t == j else 0 for t in range(k)] for j in range(k)]
    r = 0
    for c in range(m):
        while True:
            nz = [i for i in range(r, k) if M[i][c] != 0]
            if not nz: break
            p = min(nz, key=lambda i: abs(M[i][c])); M[r], M[p] = M[p], M[r]
            done = True
            for i in range(r + 1, k):
                if M[i][c] != 0:
                    f = M[i][c] // M[r][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
                    if M[i][c] != 0: done = False
            if done: r += 1; break
    return [np.array(row[m:], dtype=object) for row in M if all(v == 0 for v in row[:m])]
def lll(Bm, delta=0.99):
    """float LLL on rows of Bm (k x k), returns unimodular integer U with rows of U Bm reduced"""
    k = Bm.shape[0]; B = Bm.astype(float).copy(); U = np.eye(k, dtype=object)
    def gs():
        Q = np.zeros_like(B); mu = np.zeros((k, k))
        for i in range(k):
            Q[i] = B[i]
            for j in range(i):
                mu[i, j] = B[i] @ Q[j] / (Q[j] @ Q[j]); Q[i] -= mu[i, j] * Q[j]
        return Q, mu
    i = 1
    while i < k:
        Q, mu = gs()
        for j in range(i - 1, -1, -1):
            q = int(round(mu[i, j]))
            if q: B[i] -= q * B[j]; U[i] = U[i] - q * U[j]; Q, mu = gs()
        if Q[i] @ Q[i] >= (delta - mu[i, i - 1] ** 2) * (Q[i - 1] @ Q[i - 1]): i += 1
        else: B[[i, i - 1]] = B[[i - 1, i]]; U[[i, i - 1]] = U[[i - 1, i]]; i = max(i - 1, 1)
    return U
def short(G, bound):
    """all nonzero integer y with y^T G y <= bound, G positive definite (float)"""
    k = G.shape[0]; R = np.linalg.cholesky(G).T; out = []
    def rec(i, y, rem):
        if i < 0:
            if any(y): out.append(list(y))
            return
        c = sum(R[i, j] * y[j] for j in range(i + 1, k)); r = np.sqrt(max(rem, 0)) / R[i, i]; ctr = -c / R[i, i]
        for t in range(int(np.ceil(ctr - r - 1e-9)), int(np.floor(ctr + r + 1e-9)) + 1):
            v = (R[i, i] * t + c) ** 2
            if v <= rem + 1e-9: y[i] = t; rec(i - 1, y, rem - v)
        y[i] = 0
    rec(k - 1, [0] * k, bound); return out
def pq(S, P):
    w, V = np.linalg.eig(P)
    i = [k for k in range(N) if abs(w[k].imag) < 1e-9 and w[k].real > 1 + 1e-6][0]
    j = [k for k in range(N) if abs(w[k].imag) < 1e-9 and abs(w[k].real - 1 / w[i].real) < 1e-9][0]
    p, q = V[:, i].real, V[:, j].real; c = p @ S @ q
    return p, q / c, w[i].real
def congruent(CA, CB, slack=0.5):
    SA, SB = CA + CA.T, CB + CB.T; PA, PB = phi(CA), phi(CB)
    assert int(round(np.linalg.det(SA.astype(float)))) == int(round(np.linalg.det(SB.astype(float))))
    SAf, SBf = SA.astype(float), SB.astype(float)
    pa, qa, lam = pq(SAf, PA.astype(float)); pb, qb, lam2 = pq(SBf, PB.astype(float)); assert abs(lam - lam2) < 1e-9
    H = lambda S, p, q: S + np.outer(S @ (q - p), S @ (q - p))
    Ha, Hb = H(SAf, pa, qa), H(SBf, pb, qb); Hbi = np.linalg.inv(Hb)
    bound = N - 2 + lam + 1 / lam + slack
    # equations vec(R PB - PA R) = 0, R row-major: (I kron PB^T - PA kron I) vec(R)
    I = np.eye(N, dtype=np.int64); A = np.kron(I, PB.T) - np.kron(PA, I)
    Ks = int_kernel(A); k = len(Ks); Bs = [K.reshape(N, N) for K in Ks]; Bf = [B.astype(float) for B in Bs]
    G = np.array([[np.trace(Hbi @ Bi.T @ Ha @ Bj) for Bj in Bf] for Bi in Bf]); G = (G + G.T) / 2
    L = np.linalg.cholesky(G); U = lll(L.T); Gr = (U.astype(float) @ G @ U.astype(float).T)
    ys = short(Gr, bound); SAo, SBo = SA.astype(object), SB.astype(object); hits = []
    for y in ys:
        x = np.array(y, dtype=object) @ U     # coordinates in original kernel basis
        R = sum((x[j] * Bs[j] for j in range(k)), np.zeros((N, N), dtype=object))
        if (R.T @ SAo @ R == SBo).all(): hits.append(R)
    ok = False; CAo, CBo = CA.astype(object), CB.astype(object)
    for R in hits:
        P = R.T
        if abs(round(np.linalg.det(P.astype(float)))) == 1 and (P @ CAo @ P.T == CBo).all(): ok = True
    return dict(lam=round(lam, 4), rank=k, nvec=len(ys), nisom=len(hits), congruent=ok)
if __name__ == '__main__':
    slack = float(ARGV[1]) if len(ARGV) > 1 else 0.5
    seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
    rows = sorted(byCart.values(), key=lambda r: r['nexp']); rows = [r for r in rows if r['kept'] and not r['tp']]
    random.seed(1); CA = rows[0]['CA']; U0 = np.eye(7, dtype=np.int64)
    for _ in range(6):
        i, j = random.sample(range(7), 2); U0[i] += random.choice([-1, 1]) * U0[j]
    print('control planted:', congruent(CA, U0 @ CA @ U0.T, slack), flush=True)
    t0 = time.time(); tally = collections.Counter(); res = []
    for i, r in enumerate(rows):
        d = congruent(r['CA'], r['CB'], slack); res.append((i, r['nexp'], d)); tally[d['congruent']] += 1
        print(i, r['nexp'], d, '%.0fs' % (time.time() - t0), flush=True)
    print('TALLY', dict(tally))
