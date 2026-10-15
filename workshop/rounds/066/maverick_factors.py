"""Round 066: per-irreducible-factor Smith forms SNF(f(Phi)) (index of the Z[Phi]-module structure) for the E-190 groups at n = 11.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/066/maverick_factors.py [N]"""
import sys, collections, numpy as np, sympy
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import coxeterTables as ct, nakayama as nk
N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
x = sympy.Symbol('x')
def snf(M):
    S = smith_normal_form(sympy.Matrix(M.tolist()), domain=sympy.ZZ)
    return tuple(sorted(abs(int(S[i, i])) for i in range(S.shape[0]) if abs(int(S[i,i]))!=1))
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl in nk.allRelationLengths(N):
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
    Phi = np.rint(-Ci.T @ C).astype(np.int64); I = np.eye(N, dtype=np.int64); M = Phi.copy(); ok = False
    for m in range(200):
        if (M == I).all(): ok = True; break
        M = M @ Phi
    if ok: groups[keys[rl]].append((rl, Phi))
for key, mem in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    Phi0 = mem[0][1]; cp = sympy.Matrix(Phi0.tolist()).charpoly(x).as_expr()
    facs = sympy.factor_list(cp)[1]
    out = collections.Counter(); 
    for rl, Phi in mem:
        sig = []
        for f, e in facs:
            P = sympy.Poly(f, x); A = sympy.zeros(N, N)
            Pm = sympy.Matrix(Phi.tolist()); R = sympy.zeros(N, N); Pw = sympy.eye(N)
            for c in reversed(P.all_coeffs()): R += int(c) * Pw; Pw = Pw * Pm
            sig.append(snf(np.array(R.tolist(), dtype=object)))
        out[tuple(sig)] += 1
    print(f"group {len(mem)}: factors {[str(f) for f,e in facs]}; distinct per-factor SNF tuples {len(out)}: {dict(out)}")
