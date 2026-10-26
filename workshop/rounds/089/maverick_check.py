"""Round 089 maverick: Serre-duality check of the bimodule method: dim H^*(Hom(P_i, nu^k P_l)) = dim H^*(Hom(nu^(k-1) P_l, P_i)),
and d^2 = 0 on Hom complexes.  usage: .venv/bin/python workshop/rounds/089/maverick_check.py RELSTRING K"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from maverick_bimod import *
rl = [int(c) for c in sys.argv[1]]; K = int(sys.argv[2])
A = Alg(rl); levels = resolve(A)
cur = {}
for l in range(1, A.n + 1):
    X = Cx(); X.vert[0] = (0, l); cur[l] = X
bad = 0; tot = 0
for k in range(1, K + 1):
    prev = cur; cur = {l: nu_bimod(A, levels, X)[0] for l, X in prev.items()}
    for l in cur:
        for i in range(1, A.n + 1):
            a = hom_total(A, cur[l], i, False); b = hom_total(A, prev[l], i, True); tot += 1
            if a != b: bad += 1
    print("k", k, "R", sum(len(X.vert) for X in cur.values()), "mismatches so far", bad, "of", tot, flush=True)
