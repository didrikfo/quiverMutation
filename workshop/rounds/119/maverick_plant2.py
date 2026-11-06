"""Round 119 maverick: planted errors by construction on a low-lam n = 12 LNA Cartan matrix.  C' = C + sum of K random +-1 entry perturbations (any positions).
Keep C' in E-269's domain (det C' = 1, S' sig (11,1), det S' = det S, same Phi charpoly, one Salem pair).  Fingerprint q_S (histogram of x^T S x mod 2 on the discriminant group);
a fingerprint different from C's is a pair INEQUIVALENT BY CONSTRUCTION (Nikulin 1.9.4).  Run E-269 (maverick_iso140.test) on (C, C'): congruent False = planted error caught.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_plant2.py group(0..5 as in rounds/107/maverick_n12low.py) K tries maxtest seed"""
import sys, collections, time, random; sys.path.insert(0, "workshop/rounds/099"); sys.path.insert(0, "workshop/rounds/097")
import numpy as np, sympy as sp
from fractions import Fraction
import maverick_iso140 as M
from maverick_cartanZ import disc_form
from quivermutation import coxeterTables as ct
import os; N = int(os.environ.get("NN", 12)); M.N = N
gi, K, tries, maxtest, seed = [int(x) for x in sys.argv[1:6]]; random.seed(seed)
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
lam, mem = sel[gi]; a = mem[0]; Ca = np.array(ct.lnaCartanMatrix(N, list(a))).astype(int)
def charp(C): return tuple(int(round(x)) for x in np.poly(-np.linalg.inv(C.T.astype(float)) @ C.astype(float)))
def fp(C):
    S = sp.Matrix((C + C.T).tolist()); d, E, G, bil, q, add = disc_form(S)
    Sf = [[Fraction(int(S[i, j])) for j in range(N)] for i in range(N)]
    return tuple(sorted(collections.Counter(sum(x[i] * Sf[i][j] * x[j] for i in range(N) for j in range(N)) % 2 for x in E).items()))
cpa = charp(Ca); dSa = round(np.linalg.det(Ca + Ca.T)); fa = fp(Ca)
print("group", gi, "lam %.4f" % lam, "a", a, "detS", dSa, flush=True)
cat = collections.Counter(); found = []; seen = set(); same = []
for _ in range(tries):
    C = Ca.copy()
    for _k in range(K):
        i, j = random.randrange(N), random.randrange(N); C[i, j] += random.choice([-1, 1])
    key = C.tobytes()
    if key in seen: continue
    seen.add(key)
    if round(np.linalg.det(C)) != 1: cat["det C"] += 1; continue
    if round(np.linalg.det(C + C.T)) != dSa: cat["det S"] += 1; continue
    if (np.linalg.eigvalsh(C + C.T) < 0).sum() != 1: cat["sig"] += 1; continue
    if charp(C) != cpa: cat["charpoly"] += 1; continue
    if fp(C) == fa: cat["in domain, same q_S"] += 1; same.append(C) if len(same) < 4 else None; continue
    cat["in domain, DIFFERENT q_S"] += 1; found.append(C)
    if len(found) >= maxtest: break
print("filter counts", dict(cat), flush=True)
orig = M.cartan; cur = {}
M.cartan = lambda N_, rl: cur[rl]; cur["a"] = sp.Matrix(Ca.tolist())
for C in same:
    cur["b"] = sp.Matrix(C.tolist()); r = M.test("a", "b", 0.5)
    print("CONTROL same-q_S (expect YES):", "congruent", r["congruent"], "nsol", r["nsol"], "nvec", r["nvec"], flush=True)
for C in found:
    cur["b"] = sp.Matrix(C.tolist()); t = time.time()
    r = M.test("a", "b", 0.5); r2 = M.test("a", "b", 20)
    print("PLANTED q_S-different:", "diff entries", [(int(i), int(j), int(Ca[i, j]), int(C[i, j])) for i, j in zip(*np.nonzero(C - Ca))], "slack0.5", {k: v for k, v in r.items() if k != "P"}, "slack20 congruent", r2["congruent"], "nvec", r2["nvec"], "%.0fs" % (time.time() - t), flush=True)
