"""Round 105 maverick: scan lnaKeys groups (n=6..12) for Phi with unique real eigenvalue lam in (1,1.5), S signature (n-1,1), det C = 1;
list members and quipus; with arg 'test' run the fixed E-269 search on within-quipu pairs (a != b) and on a-vs-a.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/105/maverick_lowlam.py [scan|test]"""
import sys, collections; sys.path.insert(0, "workshop/rounds/099"); sys.path.insert(0, "workshop/rounds/097")
import numpy as np
from quivermutation import coxeterTables as ct
from quivermutation.quipuForms import quipuForAlmostSeparateLNA
def info(n, rl):
    C = np.array(ct.lnaCartanMatrix(n, list(rl))).astype(float)
    if round(np.linalg.det(C)) != 1: return None
    S = C + C.T; ev = np.linalg.eigvalsh(S)
    if not ((ev < -1e-9).sum() == 1): return None
    w = np.linalg.eigvals(-np.linalg.inv(C.T) @ C); big = [z.real for z in w if abs(z.imag) < 1e-9 and z.real > 1 + 1e-6]
    return big[0] if len(big) == 1 else None
mode = sys.argv[1] if len(sys.argv) > 1 else "scan"; cands = []
for n in range(6, 13):
    keys = ct.lnaKeys(n); groups = collections.defaultdict(list)
    for rl, kk in keys.items(): groups[kk].append(tuple(rl))
    for kk, mem in groups.items():
        lam = info(n, mem[0])
        if lam is None or lam >= 1.5: continue
        q = collections.defaultdict(list)
        for rl in mem:
            try: r = quipuForAlmostSeparateLNA(n, list(rl))
            except Exception: r = None
            if r is not None: q[r].append(rl)
        within = [v for v in q.values() if len(v) >= 2]
        print("n", n, "lam %.4f" % lam, "group", len(mem), "quipus", {i: len(v) for i, v in enumerate(q.values())}, "example", mem[:2], flush=True)
        for v in within[:1]:
            if n <= 11: cands.append((n, lam, v[0], v[1]))
if mode == "test":
    import maverick_iso140 as M
    for n, lam, a, b in cands:
        M.N = n
        for x, y, tag in ((a, b, "within"), (a, a, "self")):
            r = M.test(x, y, 0.5); print(tag, n, "lam %.4f" % lam, a, b, {k: v for k, v in r.items() if k != "P"}, flush=True)
if mode == "old":   # old l.16 threshold (real > 1.5) on the lam 1.1763 / 1.2304 within pairs: expect IndexError
    import maverick_iso140 as M
    def pq_old(S, P):
        w, V = np.linalg.eig(P); N = M.N
        i = [k for k in range(N) if abs(w[k].imag) < 1e-9 and w[k].real > 1.5][0]; return i
    for n, a in ((10, (2, 2, 2, 2, 2, 3, 0, 0)), (11, (2, 2, 2, 2, 2, 2, 3, 0, 0))):
        M.N = n; Pa = M.phi(M.cartan(n, a))
        try: pq_old(None, np.array(Pa.tolist(), dtype=float)); print(n, "old threshold ok")
        except IndexError: print(n, "old threshold 1.5: IndexError (no eigenvalue > 1.5)")
