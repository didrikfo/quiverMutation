"""Round 087 scholar: consistency check of scholar_rankgrowth.py.  The K_0 class of nu^k P_l (sum of (-1)^deg e_vertex) must equal
M^k e_l, with M the matrix of v_1 columns (nu is linear on K_0), and M must be -Phi (Phi = -C^-T C, or its transpose) from the Cartan matrix.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/087/scholar_rg_check.py RELSTRING [K]"""
import sys
import numpy as np
sys.path.insert(0, "workshop/rounds/087")
import scholar_rankgrowth as rg
from quivermutation import coxeterTables as ct
rl = [int(c) for c in sys.argv[1]]; K = int(sys.argv[2]) if len(sys.argv) > 2 else 6
A = rg.Alg(rl); n = A.n; memo = {}
cur = {l: rg.Cx() for l in range(1, n + 1)}
for l, X in cur.items(): X.vert[0] = (0, l)
def cls(X):
    v = np.zeros(n, dtype=np.int64)
    for d, x in X.vert.values(): v[x - 1] += (-1) ** d
    return v
vs = [{l: cls(X) for l, X in cur.items()}]
for k in range(1, K + 1):
    cur = {l: rg.nu_step(A, X, memo, 10**6)[0] for l, X in cur.items()}
    vs.append({l: cls(X) for l, X in cur.items()})
M = np.array([vs[1][l] for l in range(1, n + 1)]).T
ok = all((np.linalg.matrix_power(M, k) @ np.eye(n, dtype=np.int64)[:, l - 1] == vs[k][l]).all() for k in range(K + 1) for l in range(1, n + 1))
C = np.array(ct.lnaCartanMatrix(n, rl), float); Ci = np.linalg.inv(C)
Phi = np.rint(-Ci.T @ C).astype(np.int64)
print("K0 classes follow M^k:", ok, "| M == -Phi:", (M == -Phi).all(), "| M == -Phi^T:", (M == -Phi.T).all(), "| M == -Phi^-1:", (M == -np.rint(np.linalg.inv(Phi)).astype(np.int64)).all())
