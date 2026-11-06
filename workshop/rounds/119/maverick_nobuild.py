"""Round 119 maverick: build pairs INEQUIVALENT BY CONSTRUCTION inside E-269's domain, then run E-269 on them.
Random unit upper-triangular integer C (n, entries in {-1,0,1,2} sparse).  Domain: S = C + C^T of signature (n-1,1), exactly one real Phi-eigenvalue > 1 (Salem pair), lam < LAMMAX.
Collide on (charpoly of Phi, det S).  A colliding pair with q_S (x^T S x mod 2 on the discriminant group, maverick_cartanZ.disc_iso) NOT isomorphic cannot be Z-congruent
(Z-congruence of C gives isometry of S gives iso of discriminant forms, Nikulin 1.9.4) -- an independent certificate.  Then E-269 (maverick_iso140.test) must say congruent False.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_nobuild.py n tries lammax seed"""
import sys, collections, time, random; sys.path.insert(0, "workshop/rounds/099"); sys.path.insert(0, "workshop/rounds/097")
import numpy as np, sympy as sp
import maverick_iso140 as M
from maverick_cartanZ import disc_iso
n = int(sys.argv[1]); tries = int(sys.argv[2]); lammax = float(sys.argv[3]); random.seed(int(sys.argv[4])); np.random.seed(int(sys.argv[4]))
M.N = n; store = {}
def dom(C):
    S = C + C.T; d = round(np.linalg.det(S))
    if d == 0: return None
    ev = np.linalg.eigvalsh(S)
    if (ev < 0).sum() != 1: return None
    w = np.linalg.eigvals(-np.linalg.inv(C.T) @ C)
    big = [z.real for z in w if abs(z.imag) < 1e-9 and z.real > 1 + 1e-6]
    if len(big) != 1 or big[0] >= lammax: return None
    cp = tuple(int(round(x)) for x in np.poly(-np.linalg.inv(C.T) @ C))
    return (cp, abs(d)), big[0]
groups = collections.defaultdict(list)
for _ in range(tries):
    C = np.eye(n, dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.35: C[i, j] = random.choice([-1, 1, 1, 2, 0])
    r = dom(C.astype(float))
    if r: groups[r[0]].append((r[1], C))
print("domain groups", len(groups), "with >=2:", sum(len(v) >= 2 for v in groups.values()), flush=True)
orig = M.cartan; cur = {}
M.cartan = lambda N_, rl: cur[rl]
nfound = 0
for key, mem in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    seen = []
    for lam, C in mem:
        if any((C == D).all() for _, D in seen): continue
        seen.append((lam, C))
    for (l1, C1), (l2, C2) in [(seen[i], seen[j]) for i in range(len(seen)) for j in range(i + 1, len(seen))][:6]:
        A, B = sp.Matrix(C1.tolist()), sp.Matrix(C2.tolist())
        try: iso, why = disc_iso(A + A.T, B + B.T)
        except Exception as e: continue
        if iso: continue
        cur["a"], cur["b"] = A, B
        t = time.time()
        try: r = M.test("a", "b", 0.5); r2 = M.test("a", "b", 20)
        except Exception as e: print("ERR", str(e)[:80]); continue
        print("PAIR lam %.4f detS %d q_S differ(%s)" % (l1, key[1], why), "| slack0.5:", {k: v for k, v in r.items() if k != "P"}, "| slack20:", {k: v for k, v in r2.items() if k != "P"}, "| %.0fs" % (time.time() - t), flush=True)
        print(" Ca", C1.tolist(), "\n Cb", C2.tolist(), flush=True); nfound += 1
        break
    if nfound >= 4: break
