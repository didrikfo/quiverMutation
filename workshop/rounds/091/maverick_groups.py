"""Round 091 maverick: n = 11 key groups with indefinite members whose Phi has finite order or Phi^a = -I (a <= 200).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/091/maverick_groups.py [N] ; prints groups, pickles to /tmp/mav/g<N>.pkl"""
import sys, collections, pickle
import numpy as np
from quivermutation import coxeterTables as ct
N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl, k in keys.items(): groups[k].append(rl)
I = np.eye(N, dtype=np.int64); out = {}
for k, v in groups.items():
    rows = []
    for rl in v:
        C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
        neg = int((np.linalg.eigvalsh(Ci + Ci.T) < -1e-9).sum())
        if neg == 0: continue
        Phi = np.rint(-Ci.T @ C).astype(np.int64); P = Phi.copy(); o = 0; sgn = 0
        for m in range(1, 201):
            if (P == I).all(): o, sgn = m, 1; break
            if (P == -I).all(): o, sgn = m, -1; break
            P = P @ Phi
        if o: rows.append((rl, o, sgn, neg))
    if rows:
        out[k] = rows
        print("key group size", len(v), "indef periodic(+-I)", len(rows), "(order,sign,neg)", dict(collections.Counter((o, s, n) for _, o, s, n in rows)))
pickle.dump(out, open(f"/tmp/mav/g{N}.pkl", "wb"))
