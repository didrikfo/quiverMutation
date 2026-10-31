"""Round 102 skeptic response: cross-pair (a != b) test of the t-window bound.
b = P a P^T (P random unimodular), R = P^T is an exact S-isometry with R Phi_b = Phi_a R. Minimise tr(H_b^-1 R^T H_a R) over R Phi_b^k, k=-40..40
and compare with n - 2 + lam + 1/lam. Also: the 1.5 threshold in maverick_iso140.py l.16 is NOT used here (lam picked as the unique real eig > 1.0001).
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/102/skeptic_cross.py"""
import random, numpy as np
from quivermutation import coxeterTables as ct
random.seed(7); rng = np.random.default_rng(7)
def sigf(S):
    ev = np.linalg.eigvalsh(S); return int((ev > 1e-9).sum()), int((ev < -1e-9).sum())
def pqH(S, Ph):
    w, V = np.linalg.eig(Ph); n = len(w)
    re = [k for k in range(n) if abs(w[k].imag) < 1e-9 and w[k].real > 1.0001]
    if len(re) != 1: return None
    i = re[0]; lam = w[i].real
    j = [k for k in range(n) if abs(w[k].imag) < 1e-9 and abs(w[k].real - 1 / lam) < 1e-7][0]
    p, q = V[:, i].real, V[:, j].real; q = q / (p @ S @ q); v = S @ (q - p)
    return lam, S + np.outer(v, v)
def unimod(n):
    P = np.eye(n)
    for _ in range(25):
        i, j = rng.choice(n, 2, replace=False); P[i] += rng.choice([-1, 1]) * P[j]
    return P
KK = 40
tot = 0; viol = 0; worst = -9; minlam = 9; ratios = []
for n in (9, 10, 11, 12):
    keys = list(ct.lnaKeys(n).keys()); random.shuffle(keys); seen = 0
    for rl in keys:
        Ca = np.array(ct.lnaCartanMatrix(n, rl)).astype(float); Sa = Ca + Ca.T
        if sigf(Sa) != (n - 1, 1) or abs(np.linalg.det(Sa)) < 0.5: continue
        P = unimod(n); Cb = P @ Ca @ P.T; Sb = Cb + Cb.T
        Pha = -np.linalg.inv(Ca.T) @ Ca; Phb = -np.linalg.inv(Cb.T) @ Cb
        ra, rb = pqH(Sa, Pha), pqH(Sb, Phb)
        if ra is None or rb is None: continue
        lam, Ha = ra; Hbi = np.linalg.inv(rb[1]); R0 = P.T
        assert np.abs(R0.T @ Sa @ R0 - Sb).max() < 1e-6 * max(1, np.abs(R0).max()**2)
        assert np.abs(R0 @ Phb - Pha @ R0).max() < 1e-6 * max(1, np.abs(R0).max()*np.abs(Phb).max())
        best = min(np.trace(Hbi @ (R0 @ np.linalg.matrix_power(Phb, k) if k >= 0 else R0 @ np.linalg.matrix_power(np.linalg.inv(Phb), -k)).T @ Ha @
                            (R0 @ np.linalg.matrix_power(Phb, k) if k >= 0 else R0 @ np.linalg.matrix_power(np.linalg.inv(Phb), -k))) for k in range(-KK, KK + 1))
        bd = n - 2 + lam + 1 / lam; tot += 1; viol += best > bd + 1e-6 * bd; ratios.append(best - bd); minlam = min(minlam, lam)
        seen += 1
        if seen >= 40: break
print("cross pairs", tot, "violations (min over k > bound)", viol, "max(min_k value - bound)", max(ratios), "smallest lam", round(minlam, 4))
