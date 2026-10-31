"""Round 102 skeptic: re-derive E-269's search bound. Claim: BOUND - slack = n - 2 + lam + 1/lam (closed form), for ANY lam > 1.
Numpy only. For LNA Cartan matrices C (n=10..12) with S=C+C^T of signature (n-1,1): compute lam, p, q, H, the code's bound
(max over t grid of t^2.. form) and compare with closed form; check on R = Phi^k (an S-isometry commuting with Phi) that
tr(H^-1 R^T H R) = n-2 + t^2 + t^-2 with t = lam^k.
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/102/skeptic_bound.py"""
import sys, random, numpy as np
from quivermutation import coxeterTables as ct
random.seed(1)
def sigf(S):
    ev = np.linalg.eigvalsh(S); return int((ev > 1e-9).sum()), int((ev < -1e-9).sum())
tot = 0; worst = 0; worstR = 0; lams = []
for n in (10, 11, 12):
    keys = list(ct.lnaKeys(n).keys()); random.shuffle(keys); seen = 0
    for rl in keys:
        C = np.array(ct.lnaCartanMatrix(n, rl)).astype(float); S = C + C.T
        if sigf(S) != (n - 1, 1) or abs(np.linalg.det(S)) < 0.5: continue
        Ph = -np.linalg.inv(C.T) @ C
        w, V = np.linalg.eig(Ph)
        re = [k for k in range(n) if abs(w[k].imag) < 1e-9 and w[k].real > 1.0001]
        if len(re) != 1: continue
        i = re[0]; lam = w[i].real
        j = [k for k in range(n) if abs(w[k].imag) < 1e-9 and abs(w[k].real - 1 / lam) < 1e-7][0]
        p, q = V[:, i].real, V[:, j].real; q = q / (p @ S @ q)
        v = S @ (q - p); H = S + np.outer(v, v); Hi = np.linalg.inv(H)
        assert np.linalg.eigvalsh(H).min() > 0
        ts = np.exp(np.linspace(-0.5 * np.log(lam), 0.5 * np.log(lam), 2001))
        code = np.trace(Hi @ S) + max((t * q - p / t) @ S @ Hi @ S @ (t * q - p / t) for t in ts)
        closed = n - 2 + lam + 1 / lam
        worst = max(worst, abs(code - closed))
        for k in (0, 1, -1, 2):
            R = np.linalg.matrix_power(Ph, k) if k >= 0 else np.linalg.matrix_power(np.linalg.inv(Ph), -k)
            assert np.abs(R.T @ S @ R - S).max() < 1e-6 * max(1, np.abs(R).max() ** 2)
            lhs = np.trace(Hi @ R.T @ H @ R); rhs = n - 2 + lam ** (2 * k) + lam ** (-2 * k)
            worstR = max(worstR, abs(lhs - rhs) / rhs)
        lams.append((round(lam, 4), n)); tot += 1; seen += 1
        if seen >= 150: break
print("matrices checked", tot, "max |code bound - closed form|", worst, "max rel err of identity on Phi^k", worstR)
lo = sorted(set(lams))[:6]; print("smallest lam (lam, n):", lo)
for lam in (1.1, 1.17628, 1.3, 1.5, 2.1451):
    print("lam", lam, "bound without slack n=12:", 10 + lam + 1 / lam)
