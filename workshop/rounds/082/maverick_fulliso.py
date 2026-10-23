"""Round 082 (idea): FULL-lattice exact equivariant isometry for the INDEFINITE n = 11 15-group (stage 2+3 of r074 in one go, no number theory).
Idea: G = C^-1 + C^-T is M-invariant (M = Phi^T), signature (9,2); the 2 negative directions are the M-eigenplane of angle 1/14 (E-209 sector Phi_14).
H = G (I - 2 pi_N) is an M-invariant POSITIVE DEFINITE majorant; an M-equivariant isometry P: G_a -> G_b maps N_a -> N_b (same eigenvalue, checked by
maverick_eigsig.py: one profile for all 15), hence is also an H-isometry, so P v_i lies in the FINITE set {u : H_b(u) = H_a(v_i)} (Fincke-Pohst).
Same generator/orbit scheme as r074 maverick_isometry.isometry; the final P is checked EXACTLY (integer, det +-1, P^T G_b P = G_a, P M_a = M_b P).
A NO is a proof of non-isometry (hence Cartan non-congruence) only if the enumeration is exhaustive: float tolerance 1e-6, node cap reported.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/082/maverick_fulliso.py [--plan] [--n10]"""
import sys, time, collections
sys.path.insert(0, "workshop/rounds/074"); sys.path.insert(0, "workshop/rounds/082")
import numpy as np
from sympy import Matrix
from maverick_isometry import data, shortVectors, orbitCols, cols_for
from maverick_eigsig import isoClasses15

def majorant(G, M, negAngle):
    """H = G (I - 2 pi_N), N = sum of M-eigenspaces with eigenvalue angle in negAngle (set of angles/2pi, rounded 4)."""
    lam, W = np.linalg.eig(M.astype(float)); n = len(lam); Gf = G.astype(float); pi = np.zeros((n, n), complex)
    for k in range(n):
        if round((np.angle(lam[k]) / (2 * np.pi)) % 1.0, 4) in negAngle:
            w = W[:, k:k+1]; pi += (w @ w.conj().T @ Gf) / (w.conj().T @ Gf @ w)[0, 0]
    assert abs(pi.imag).max() < 1e-8; pi = pi.real
    H = Gf @ (np.eye(n) - 2 * pi); assert np.allclose(H, H.T, atol=1e-8)
    assert np.linalg.eigvalsh((H + H.T) / 2).min() > 1e-6, "majorant not positive definite"
    return (H + H.T) / 2

def isometryH(Ga, Ma, Gb, Mb, Ha, Hb, cap=2000000, kmax=5, tol=1e-6):
    n = Ga.shape[0]
    Mp = [np.linalg.matrix_power(Ma, c) for c in range(n)]; Mq = [np.linalg.matrix_power(Mb, c) for c in range(n)]
    def span(ws): return np.linalg.matrix_rank(np.array([Mp[c] @ w for w in ws for c in range(n)], dtype=float)) if ws else 0
    hn = lambda H, v: float(v @ H @ v)
    vs = []
    for v in sorted([v for v in shortVectors(Ha, 3.0) if v.any()], key=lambda z: hn(Ha, z)):
        if span(vs + [v]) > span(vs): vs = vs + [v]
        if span(vs) == n: break
    if span(vs) < n or len(vs) > kmax: return "SKIP-gens", None, 0, vs
    idx, W = orbitCols(Ma, vs, n); Winv = W.inv()
    Bmax = max(hn(Ha, v) for v in vs)
    byNorm = collections.defaultdict(list)
    for u in shortVectors(Hb, Bmax):
        if u.any(): byNorm[int(u @ Gb @ u)].append(u)
    nodes = [0]; found = [None]
    def rec(k, us):
        if found[0] is not None or nodes[0] > cap: return
        if k == len(vs):
            P = cols_for(Mb, us, idx) * Winv
            if any(x.q != 1 for x in P) or abs(P.det()) != 1: return
            if P.T * Matrix(Gb.tolist()) * P == Matrix(Ga.tolist()) and P * Matrix(Ma.tolist()) == Matrix(Mb.tolist()) * P: found[0] = P
            return
        q = int(vs[k] @ Ga @ vs[k])
        for u in byNorm[q]:
            if abs(hn(Hb, u) - hn(Ha, vs[k])) > tol: continue
            nodes[0] += 1
            ok = all(int(vs[k] @ Ga @ Mp[c] @ vs[j]) == int(u @ Gb @ Mq[c] @ us[j]) and int(vs[j] @ Ga @ Mp[c] @ vs[k]) == int(us[j] @ Gb @ Mq[c] @ u)
                     for j in range(k) for c in range(n)) and all(int(vs[k] @ Ga @ Mp[c] @ vs[k]) == int(u @ Gb @ Mq[c] @ u) for c in range(n))
            if ok: rec(k + 1, us + [u])
    rec(0, [])
    if found[0] is not None: return "YES", found[0], nodes[0], vs
    return ("INCONCLUSIVE-cap" if nodes[0] > cap else "NO"), None, nodes[0], vs

def run(pairs, N, neg, label, t0):
    for a, b in pairs:
        Ga, Ma, _ = data(a, N); Gb, Mb, _ = data(b, N)
        Ha, Hb = majorant(Ga, Ma, neg), majorant(Gb, Mb, neg)
        res, P, nodes, vs = isometryH(Ga, Ma, Gb, Mb, Ha, Hb)
        Ca = np.array(__import__("quivermutation.coxeterTables", fromlist=["x"]).lnaCartanMatrix(N, tuple(a)), dtype=np.int64)
        Cb = np.array(__import__("quivermutation.coxeterTables", fromlist=["x"]).lnaCartanMatrix(N, tuple(b)), dtype=np.int64)
        cart = ""
        if P is not None:
            Pn = np.array(P.tolist(), dtype=np.int64); cart = " Cartan P C_a P^T == C_b: %s" % bool((Pn @ Ca @ Pn.T == Cb).all())
        print(label, "".join(map(str, a)), "->", "".join(map(str, b)), res, "gens", len(vs), "nodes", nodes, cart, f"{time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    t0 = time.time(); rd = lambda s: tuple(int(c) for c in s)
    if "--plan" in sys.argv:
        print("plan: n=10 power controls (50505000 vs 34504030 indefinite, YES expected? E-178; D_10 22222230 vs 50505000 should be NO by signature -> majorant undefined, skipped) ~seconds;")
        print("n=11: 8 pairs rep0 vs rep_i, each Fincke-Pohst on a rank-11 majorant to H-norm <= ~3: sizes printed; load derivedOrbits ~35 s"); sys.exit()
    # n = 10 controls: equivalent indefinite pair 50505000 ~ 34504030?  (listed as same class in E-190: check by running); self pair is a sanity check
    neg10 = {0.0556, 0.9444}
    run([(rd("50505000"), rd("50505000")), (rd("50505000"), rd("34504030")), (rd("34504030"), rd("50505000"))], 10, neg10, "n10", t0)
    neg11 = {0.0714, 0.9286}
    mem, reps, find = isoClasses15()
    run([(reps[0], r) for r in reps], 11, neg11, "n11", t0)
