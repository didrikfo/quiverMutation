"""Round 110 skeptic: for n = 8..12, every lnaKeys group with >= 2 quipus (first 400 members each, almost-separate members only): lam of Phi, signature, #quipus.
Question: does a q_S-separable (>= 2 quipu) group exist with lam < 1.5?  usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/110/skeptic_lamscan.py"""
import sys, collections; sys.path.insert(0, "workshop/rounds/097")
import numpy as np
from maverick_cartanZ import *
for N in range(8, 13):
    keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
    for rl, k in keys.items(): groups[k].append(rl)
    lams = []; ng = 0
    for k, v in groups.items():
        qs = {quipuForAlmostSeparateLNA(N, list(rl)) for rl in v}; qs.discard(None)
        ng += 1
        if len(qs) >= 2:
            C = cartan(N, v[0]); lam = max(abs(np.linalg.eigvals(np.array(phi(C).tolist(), dtype=float)))); lams.append(round(lam, 4))
    print(N, "groups", ng, "with>=2 quipus", len(lams), "min lam", min(lams) if lams else None, sorted(lams)[:5], flush=True)
