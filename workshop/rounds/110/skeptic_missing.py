"""Round 110 skeptic: which >=2-quipu n=12 group does skeptic_qsno.py skip, and why (domain filter)."""
import sys, collections; sys.path.insert(0, "workshop/rounds/097")
import numpy as np
from maverick_cartanZ import *
N = 12; keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl, k in keys.items(): groups[k].append(rl)
for k, v in groups.items():
    q = collections.defaultdict(list)
    for rl in v:
        r = quipuForAlmostSeparateLNA(N, list(rl))
        if r is not None: q[r].append(tuple(rl))
    if len(q) < 2: continue
    qs = list(q); a = q[qs[0]][0]; b = q[qs[1]][0]; Ca, Cb = cartan(N, a), cartan(N, b); Sa, Sb = Ca + Ca.T, Cb + Cb.T
    lam = max(abs(np.linalg.eigvals(np.array(phi(Ca).tolist(), dtype=float))))
    print(len(v), round(lam, 4), sig(Sa), sig(Sb), int(Sa.det()), int(Sb.det()), int(Ca.det()))
