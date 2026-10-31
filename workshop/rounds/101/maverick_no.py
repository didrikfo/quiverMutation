"""Round 101 maverick: lowest-lam NO control available.  Census n=10,11,12 (maverick_find.py): every same-key class with lam<1.5 has ONE quipu, so no certified-inequivalent
pair with lam<1.5 exists; the smallest cross-quipu lam is 1.66 (n=12, det -12).  Run E-269 test (patched pq, threshold 1.05) on cross-quipu pairs of that class; also a
synthetic small-lam NO: N(10,8) vs a form S' with the same Phi-class but S' = S_b of a non-congruent twin is not available, so only the lam=1.66 pair.
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/101/maverick_no.py"""
import sys, collections; sys.path.insert(0, "workshop/rounds/099"); sys.path.insert(0, "workshop/rounds/097")
import maverick_iso140 as M
from maverick_cartanZ import *
import numpy as np
def pq(S, P):
    N = M.N; w, V = np.linalg.eig(P)
    i = [k for k in range(N) if abs(w[k].imag) < 1e-9 and w[k].real > 1.05][0]; j = [k for k in range(N) if abs(w[k].imag) < 1e-9 and abs(w[k].real - 1 / w[i].real) < 1e-9][0]
    p, q = V[:, i].real, V[:, j].real; c = p @ S @ q
    return p, q / c, w[i].real
M.pq = pq; n = 12; M.N = n
keys = ct.lnaKeys(n); groups = collections.defaultdict(list)
for rl, kk in keys.items(): groups[kk].append(rl)
for kk, g in groups.items():
    C = cartan(n, g[0]); S = C + C.T
    if sig(S) != (n - 1, 1) or int(S.det()) != -12: continue
    q = collections.defaultdict(list)
    for rl in g:
        r = quipuForAlmostSeparateLNA(n, list(rl))
        if r is not None: q[r].append(rl)
    if len(q) < 2: continue
    reps = [v[0] for v in q.values()]; print("quipus", list(q), "reps", reps, flush=True)
    for slack in (0.5, 30):
        r = M.test(reps[0], reps[1], slack); print("slack", slack, {k: v for k, v in r.items() if k != "P"}, flush=True)
    r = M.test(q[list(q)[0]][0], q[list(q)[0]][1], 0.5); print("within-quipu control", {k: v for k, v in r.items() if k != "P"}, flush=True)
