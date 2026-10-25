"""Round 087 scholar: list indefinite LNAs at n with Phi^a = +-I (period a) or spectral radius > 2.5.
usage (repo root): .venv/bin/python workshop/rounds/087/scholar_rg_list.py N     -> lines: relstring period rho"""
import sys, numpy as np
from quivermutation import coxeterTables as ct, nakayama as nk
N = int(sys.argv[1])
for rl in nk.allRelationLengths(N):
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
    Phi = np.rint(-Ci.T @ C).astype(np.int64); I = np.eye(N, dtype=np.int64); M = Phi.copy(); o = 0
    for m in range(1, 300):
        if (M == I).all() or (M == -I).all(): o = m; break
        M = M @ Phi
    rho = np.abs(np.linalg.eigvals(Phi)).max()
    if o or rho > (float(sys.argv[2]) if len(sys.argv) > 2 else 2.5): print("".join(map(str, rl)), o, round(rho, 3))
