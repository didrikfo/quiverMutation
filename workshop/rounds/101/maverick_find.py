"""Round 101 maverick: find cross-quipu (certified inequivalent) LNA pairs at n=10,11 with equal Phi charpoly, sig S=(n-1,1), det S !=0, Salem lam<1.5.
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/101/maverick_find.py"""
import sys, collections, itertools; sys.path.insert(0, "workshop/rounds/097")
from maverick_cartanZ import *
for n in (10, 11):
    keys = ct.lnaKeys(n); groups = collections.defaultdict(list)
    for rl, kk in keys.items(): groups[kk].append(rl)
    print(n, "groups", len(groups), "strings", len(keys), flush=True)
    for kk, g in groups.items():
        if len(g) < 2: continue
        C = cartan(n, g[0]); S = C + C.T
        if sig(S) != (n - 1, 1) or int(S.det()) == 0: continue
        w = np.linalg.eigvals(np.array(phi(C).tolist(), dtype=float)); lam = max(w.real[abs(w.imag) < 1e-9])
        q = collections.defaultdict(list)
        for rl in g:
            r = quipuForAlmostSeparateLNA(n, list(rl))
            if r is not None: q[r].append(rl)
        print(n, "group", len(g), "det", int(S.det()), "lam %.4f" % lam, "quipus", {str(k): len(v) for k, v in q.items()}, flush=True)
