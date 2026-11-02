"""Round 109 maverick: planted-error control, two-entry perturbations.  Single-entry in-domain perturbations of a low-lam C were all Z-congruent to C
(see maverick_planted.py), so they are not errors.  Here C' = C + d1 E_ij + d2 E_kl (d = +-1) is filtered to the test domain (det C' = 1, S' signature (11,1),
det S' = det S, same Phi charpoly, found in floating point then exact) and E-269 (maverick_iso140.test) is run on (C, C') for up to MAXT candidates.
Every candidate with congruent False is a planted error caught (NO).  usage (repo root):
timeout 10m .venv/bin/python -W ignore workshop/rounds/109/maverick_planted2.py [group = 5] [MAXT = 12] [slack = 0.5]"""
import sys, collections, time, random; sys.path.insert(0, "workshop/rounds/099"); sys.path.insert(0, "workshop/rounds/097")
import numpy as np, sympy as sp
import maverick_iso140 as M
from quivermutation import coxeterTables as ct
N = 12; M.N = N
gi = int(sys.argv[1]) if len(sys.argv) > 1 else 5; MAXT = int(sys.argv[2]) if len(sys.argv) > 2 else 12; slack = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
def lam_of(rl):
    C = np.array(ct.lnaCartanMatrix(N, list(rl))).astype(float)
    if round(np.linalg.det(C)) != 1: return None
    ev = np.linalg.eigvalsh(C + C.T)
    if (ev < -1e-9).sum() != 1: return None
    w = np.linalg.eigvals(-np.linalg.inv(C.T) @ C); big = [z.real for z in w if abs(z.imag) < 1e-9 and z.real > 1 + 1e-6]
    return big[0] if len(big) == 1 else None
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl, kk in keys.items(): groups[kk].append(tuple(rl))
sel = sorted([(lam_of(m[0]), m) for m in groups.values() if lam_of(m[0]) and lam_of(m[0]) < 1.5], key=lambda t: t[0])
lam, mem = sel[gi]; a = mem[0]; print("group", gi, "lam %.6f" % lam, "a", a, flush=True)
Cn = np.array(ct.lnaCartanMatrix(N, list(a))).astype(int); Ca = sp.Matrix(Cn.tolist()); orig = M.cartan
Cb = None
M.cartan = lambda n, rl: Cb if rl == "PERT" else orig(n, rl)
def phip(C): return np.poly(-np.linalg.inv(C.T.astype(float)) @ C.astype(float)).real
pa = phip(Cn); detSa = round(np.linalg.det((Cn + Cn.T).astype(float)))
cells = [(i, j) for i in range(N) for j in range(N)]; random.seed(3); cand = []; tot = 0
for x in range(len(cells)):
    for y in range(x + 1, len(cells)):
        for d1 in (1, -1):
            for d2 in (1, -1):
                C = Cn.copy(); C[cells[x]] += d1; C[cells[y]] += d2; tot += 1
                if round(np.linalg.det(C.astype(float))) != 1: continue
                S = (C + C.T).astype(float)
                if round(np.linalg.det(S)) != detSa or (np.linalg.eigvalsh(S) < -1e-9).sum() != 1: continue
                if np.abs(phip(C) - pa).max() > 1e-6: continue
                cand.append((cells[x], d1, cells[y], d2))
print("two-entry perturbations", tot, "in domain", len(cand), flush=True)
random.shuffle(cand); res = collections.Counter()
for (c1, d1, c2, d2) in cand[:MAXT]:
    Cb = Ca.copy(); Cb[c1] += d1; Cb[c2] += d2; t = time.time()
    try: r = M.test(a, "PERT", slack); res["congruent=%s" % r["congruent"]] += 1; print((c1, d1, c2, d2), "nvec", r["nvec"], "nsol", r["nsol"], "congruent", r["congruent"], "%.0fs" % (time.time() - t), flush=True)
    except Exception as e: res["ERROR " + type(e).__name__] += 1; print((c1, d1, c2, d2), "ERROR", type(e).__name__, str(e)[:80], flush=True)
print("summary", dict(res))
