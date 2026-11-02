"""Round 109 maverick: planted-error control for the E-269 test (maverick_iso140.test) at n = 12, lam < 1.5.
For one low-lam lnaKeys group (index in the lam-sorted list of maverick_n12low.py) take a = first member, b = same member with ONE Cartan entry
C[i,j] += d (d = +-1, all i,j).  Each perturbed C' is classified: (A) leaves the test domain (det C' != 1, S' not sig (11,1), det S' != det S_a, or
Phi' has a different charpoly/lam) -- trivially NO, no power shown; (B) stays in the domain (same Phi charpoly, same det S) -- then E-269 must return
congruent False, and congruent True would be a bug.  Control: unperturbed self pair must be YES.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/109/maverick_planted.py [group index = 5] [slack = 0.5]"""
import sys, collections, time; sys.path.insert(0, "workshop/rounds/099"); sys.path.insert(0, "workshop/rounds/097")
import numpy as np, sympy as sp
import maverick_iso140 as M
from quivermutation import coxeterTables as ct
N = 12; M.N = N
gi = int(sys.argv[1]) if len(sys.argv) > 1 else 5; slack = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
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
lam, mem = sel[gi]; a = mem[0]; print("group", gi, "lam %.6f" % lam, "size", len(mem), "a", a, flush=True)
Ca = sp.Matrix(np.array(ct.lnaCartanMatrix(N, list(a))).astype(int).tolist()); orig = M.cartan
def patched(n, rl):
    return Cb if rl == "PERT" else orig(n, rl)
M.cartan = patched
def phichar(C): return sp.Poly((-(C.T.inv()) * C).charpoly().as_expr(), sp.Symbol("lambda"))
cpa = phichar(Ca).all_coeffs(); detSa = int((Ca + Ca.T).det())
cats = collections.Counter(); inb = []
for i in range(N):
    for j in range(N):
        for d in (1, -1):
            Cb = Ca.copy(); Cb[i, j] += d
            if Cb.det() != 1: cats["A:det C != 1"] += 1; continue
            S = Cb + Cb.T
            if int(S.det()) != detSa: cats["A:det S differs"] += 1; continue
            if phichar(Cb).all_coeffs() != cpa: cats["A:Phi charpoly differs"] += 1; continue
            inb.append((i, j, d))
print("domain filter:", dict(cats), "| in-domain perturbations (B):", inb, flush=True)
t = time.time(); r = M.test(a, a, slack); print("control self:", {k: v for k, v in r.items() if k != "P"}, "%.0fs" % (time.time() - t), flush=True)
for (i, j, d) in inb:
    Cb = Ca.copy(); Cb[i, j] += d; t = time.time()
    try: r = M.test(a, "PERT", slack); print("B", (i, j, d), {k: v for k, v in r.items() if k != "P"}, "%.0fs" % (time.time() - t), flush=True)
    except Exception as e: print("B", (i, j, d), "ERROR", type(e).__name__, str(e)[:100], flush=True)
