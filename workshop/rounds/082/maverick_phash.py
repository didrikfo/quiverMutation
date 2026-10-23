"""Round 082 response: dump reps, H-norm bounds, node caps and sha256 of every integer P (rep0 -> rep_i, n=11) to workshop/rounds/082/maverick_P.txt.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/082/maverick_phash.py"""
import sys, hashlib
sys.path.insert(0, "workshop/rounds/074"); sys.path.insert(0, "workshop/rounds/082")
import numpy as np
from maverick_fulliso import majorant, isometryH
from maverick_isometry import data
from maverick_eigsig import isoClasses15
neg = {0.0714, 0.9286}; out = []
mem, reps, find = isoClasses15()
Ga, Ma, _ = data(reps[0], 11); Ha = majorant(Ga, Ma, neg)
allP = []
for i, b in enumerate(reps):
    Gb, Mb, _ = data(b, 11); Hb = majorant(Gb, Mb, neg)
    res, P, nodes, vs = isometryH(Ga, Ma, Gb, Mb, Ha, Hb)
    Bmax = max(float(v @ Ha @ v) for v in vs)
    s = str(P.tolist()) if P is not None else "none"
    h = hashlib.sha256(s.encode()).hexdigest()
    allP.append(s)
    out.append(f"rep{i} key {''.join(map(str,b))} {res} gens {len(vs)} Hnorm_bound(Bmax) {Bmax:.4f} nodes {nodes} cap 2000000 cap_hit {nodes>2000000} sha256(P) {h}")
    print(out[-1], flush=True)
    if i == 0 or i == 1: out.append(f"P_{i} = {s}")
out.append("sha256(all nine P, concatenated repr) " + hashlib.sha256("|".join(allP).encode()).hexdigest())
out.append("rep0 = " + "".join(map(str, reps[0])))
open("workshop/rounds/082/maverick_P.txt", "w").write("\n".join(out) + "\n")
print(out[-2])
