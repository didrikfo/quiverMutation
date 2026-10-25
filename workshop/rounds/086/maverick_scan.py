"""Round 086 (idea): size/scan key-equal groups of LNAs at n = 12, 13 for a Cartan-incongruent pair (real NO power).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/086/maverick_scan.py N [--plan]
Stage A: group all LNAs of length N by Coxeter key; keep groups with >= 2 members; classify each member as definite / indefinite
(G = C^-1 + C^-T) and periodic (Phi^m = I, m <= 200 or Phi^m = -I handled by the same loop on Phi^2 order)."""
import sys, time, collections
import numpy as np
from quivermutation import coxeterTables as ct
import os, tempfile
N = int(sys.argv[1]); t0 = time.time()
keys = ct.lnaKeys(N); print("n", N, "LNAs", len(keys), f"{time.time()-t0:.0f}s", flush=True)
groups = collections.defaultdict(list)
for rl, k in keys.items(): groups[k].append(rl)
multi = {k: v for k, v in groups.items() if len(v) > 1}
print("key groups", len(groups), "with >=2 members", len(multi), "LNAs in them", sum(map(len, multi.values())), "max size", max(map(len, multi.values())))
if "--plan" in sys.argv: sys.exit()
I = np.eye(N, dtype=np.int64); stat = collections.Counter(); per = {}
for k, v in multi.items():
    for rl in v:
        C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
        ev = np.linalg.eigvalsh(Ci + Ci.T); neg = int((ev < -1e-9).sum())
        Phi = np.rint(-Ci.T @ C).astype(np.int64); P = Phi.copy(); per_ = 0
        for m in range(1, 201):
            if (P == I).all(): per_ = m; break
            P = P @ Phi
        per[rl] = (neg, per_)
cnt = collections.Counter((neg > 0, p > 0) for neg, p in per.values()); print("(indefinite, periodic) counts", dict(cnt))
sel = [(k, v) for k, v in multi.items() if any(per[r][1] > 0 and per[r][0] > 0 for r in v)]
print("groups with an indefinite periodic member", len(sel), "sizes", sorted(len(v) for k, v in sel)[-10:])
import pickle; pickle.dump({"multi": multi, "per": per}, open(os.environ.get("MAVERICK_PKL", os.path.join(tempfile.gettempdir(), "maverick_scan%d.pkl" % N)), "wb"))
print(f"{time.time()-t0:.0f}s")
