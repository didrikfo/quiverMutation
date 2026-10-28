"""Round 094 maverick: for each n = 11 (and 9, 10) key group with >= 2 quipu classes, Cartan-level probes on 2 members per quipu:
det C, signature of C+C^T, SNF(xC+C^T) for x = -5..5 (as in E-209), SNF of C+C^T, and the F-047 profile equality from maverick_quipuprofile.py.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/094/maverick_quipucartan.py N"""
import sys, collections
import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import coxeterTables as ct
from quivermutation.quipuForms import quipuForAlmostSeparateLNA
N = int(sys.argv[1]); keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl, k in keys.items(): groups[k].append(rl)
def probe(rl):
    C = sp.Matrix(np.array(ct.lnaCartanMatrix(N, rl)).astype(int).tolist())
    G = np.array((C + C.T).tolist(), float); ev = np.linalg.eigvalsh(G)
    snf = lambda M: tuple(sorted(abs(int(smith_normal_form(M, domain=sp.ZZ)[i, i])) for i in range(N)))
    return (int(C.det()), (int((ev>1e-9).sum()), int((ev<-1e-9).sum())), tuple(snf(x*C + C.T) for x in range(-5, 6)))
for k, v in groups.items():
    q = collections.defaultdict(list)
    for rl in v:
        r = quipuForAlmostSeparateLNA(N, list(rl))
        if r is not None: q[r].append(rl)
    if len(q) < 2: continue
    P = {r: [probe(rl) for rl in m[:2]] for r, m in q.items()}
    same = len({tuple(sorted(set(p))) for p in P.values()}) == 1 and all(len(set(p)) == 1 for p in P.values())
    print("size", len(v), list(q), "det", [p[0][0] for p in P.values()], "sig", [p[0][1] for p in P.values()], "all probes equal across quipus:", same)
    print("   reps", [m[0] for m in q.values()])
