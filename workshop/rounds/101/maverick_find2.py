"""Round 101 maverick: look for pairs at n=10,11 with equal Phi charpoly, equal det S, sig (n-1,1), Salem lam<1.5 but in DIFFERENT lnaKeys classes (candidate small-lam NO controls).
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/101/maverick_find2.py"""
import sys, collections; sys.path.insert(0, "workshop/rounds/097")
from maverick_cartanZ import *
x = sp.symbols("x")
for n in (10, 11):
    keys = ct.lnaKeys(n); b = collections.defaultdict(lambda: collections.defaultdict(list)); cache = {}
    for rl, kk in keys.items():
        if kk not in cache:
            C = cartan(n, rl); S = C + C.T
            if sig(S) != (n - 1, 1) or int(S.det()) == 0: cache[kk] = None
            else:
                w = np.linalg.eigvals(np.array(phi(C).tolist(), dtype=float)); lam = max(w.real[abs(w.imag) < 1e-9])
                cache[kk] = (str(phi(C).charpoly(x).as_expr()), int(S.det()), lam) if lam < 1.5 else None
        if cache[kk]: b[cache[kk][:2]][kk].append(rl)
    for k, d in b.items():
        print(n, "det", k[1], "classes", len(d), [len(v) for v in d.values()], flush=True)
        if len(d) > 1:
            for kk, v in d.items(): print("   ", v[0])
