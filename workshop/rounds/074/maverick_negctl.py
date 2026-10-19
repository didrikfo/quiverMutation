"""Round 074 NO-control for maverick_isometry: compare (G, M) of 90000000 (n = 10) with (G, M^k), k coprime to 18 (Galois twist; G is M^k-invariant).
Also vs (G', M) with G' the Gram of a different definite lattice of equal det (A_10-key 00000000 vs D_10: det G 11 vs 4 -> not equal, shown only for the log).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/074/maverick_negctl.py"""
import sys, time, math
sys.path.insert(0, "workshop/rounds/074")
import numpy as np
from maverick_isometry import data, isometry
t0 = time.time()
for s in ("90000000", "22222230"):
    rl = tuple(int(c) for c in s); G, M, _ = data(rl, 10)
    for k in (1, 5, 7, 11, 13, 17):
        Mk = np.linalg.matrix_power(M, k)
        r, P, nodes = isometry(G, M, G, Mk)
        print(s, "k =", k, r, "nodes", nodes, f"{time.time()-t0:.1f}s", flush=True)
