"""Round 094 maverick: do periodic key groups (Phi^a = +-I, a <= 200) contain two LNAs named by DIFFERENT quipus
(Ladkani-type quipu theorem, certified derived inequivalence)?  Also counts, over ALL key groups, those with >= 2 quipus.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/094/maverick_quipupairs.py N"""
import sys, collections
import numpy as np
from quivermutation import coxeterTables as ct
from quivermutation.quipuForms import quipuForAlmostSeparateLNA
N = int(sys.argv[1]); keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl, k in keys.items(): groups[k].append(rl)
I = np.eye(N, dtype=np.int64)
def periodic(rl):
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Phi = np.rint(-np.linalg.inv(C).T @ C).astype(np.int64); P = Phi.copy()
    for m in range(1, 201):
        if (P == I).all(): return m, 1
        if (P == -I).all(): return m, -1
        P = P @ Phi
    return None
allMulti = 0; perGroups = 0; perMulti = []; namedInPer = 0; perMembers = 0
for k, v in groups.items():
    q = {}
    for rl in v:
        r = quipuForAlmostSeparateLNA(N, list(rl)); 
        if r is not None: q[rl] = r
    multi = len(set(q.values())) >= 2
    allMulti += multi
    per = any(periodic(rl) for rl in v) if len(v) < 5000 else periodic(v[0])   # any member periodic
    if per:
        perGroups += 1; perMembers += len(v); namedInPer += len(q)
        if multi: perMulti.append((k, len(v), sorted(set(q.values()))))
print("N", N, "key groups", len(groups), "groups with >=2 distinct quipus", allMulti)
print("periodic (+-I, a<=200) groups", perGroups, "members", perMembers, "quipu-named members in them", namedInPer, "periodic groups with >=2 quipus", len(perMulti))
for x in perMulti: print(x)
