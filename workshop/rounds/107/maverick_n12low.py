"""Round 107 maverick: E-269 Z-congruence search (E-276 closed-form bound) on the n = 12 low-lam (lam < 1.5) lnaKeys groups.
For each group: first two same-quipu members (YES expected; pre-registered), a-vs-a self control, and for groups with no quipu
the first pairs of members (verdict only; equivalence uncertified).  usage (repo root):
timeout 10m .venv/bin/python -W ignore workshop/rounds/107/maverick_n12low.py [slack] [group index list]"""
import sys, collections, time; sys.path.insert(0, "workshop/rounds/099"); sys.path.insert(0, "workshop/rounds/097")
import numpy as np
import maverick_iso140 as M
from quivermutation import coxeterTables as ct
from quivermutation.quipuForms import quipuForAlmostSeparateLNA
N = 12; M.N = N
slack = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
only = [int(x) for x in sys.argv[2:]]
def info(rl):
    C = np.array(ct.lnaCartanMatrix(N, list(rl))).astype(float)
    if round(np.linalg.det(C)) != 1: return None
    ev = np.linalg.eigvalsh(C + C.T)
    if (ev < -1e-9).sum() != 1: return None
    w = np.linalg.eigvals(-np.linalg.inv(C.T) @ C); big = [z.real for z in w if abs(z.imag) < 1e-9 and z.real > 1 + 1e-6]
    return big[0] if len(big) == 1 else None
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl, kk in keys.items(): groups[kk].append(tuple(rl))
sel = []
for kk, mem in groups.items():
    lam = info(mem[0])
    if lam is not None and lam < 1.5: sel.append((lam, mem))
sel.sort(key=lambda t: t[0])
for gi, (lam, mem) in enumerate(sel):
    if only and gi not in only: continue
    q = collections.defaultdict(list)
    for rl in mem:
        try: r = quipuForAlmostSeparateLNA(N, list(rl))
        except Exception: r = None
        if r is not None: q[r].append(rl)
    w = np.linalg.eigvals(-np.linalg.inv(np.array(ct.lnaCartanMatrix(N, list(mem[0]))).astype(float).T) @ np.array(ct.lnaCartanMatrix(N, list(mem[0]))).astype(float))
    print("GROUP", gi, "lam %.6f" % lam, "size", len(mem), "quipus", [len(v) for v in q.values()], "ex", mem[0], flush=True)
    pairs = []
    for v in q.values():
        if len(v) >= 2: pairs.append(("within-quipu", v[0], v[1])); break
    if not q: pairs += [("noquipu", mem[0], mem[1]), ("noquipu", mem[0], mem[-1])]
    pairs.append(("self", mem[0], mem[0]))
    for tag, a, b in pairs:
        t = time.time()
        try: r = M.test(a, b, slack); print(" ", tag, a, b, {k: v for k, v in r.items() if k != "P"}, "%.0fs" % (time.time() - t), flush=True)
        except Exception as e: print(" ", tag, a, b, "ERROR", type(e).__name__, str(e)[:100], flush=True)
