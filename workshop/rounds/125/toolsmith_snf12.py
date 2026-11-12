"""Round 125 (toolsmith): sizing T10 at n = 12.  Coxeter-key groups of LNAs of length N whose members split by (SNF(C + C^T), signature);
q_S is NOT computed (it costs most of round 122's group table), so a group split only by q_S is missed.  Prints the split groups, their class
sizes, SNF tails, and a lexicographically first root per class.  usage (repo root):
  timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_snf12.py [N=12]"""
import sys, time, collections
sys.path.insert(0, '.')
import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import coxeterTables as ct
N = int(sys.argv[1]) if len(sys.argv) > 1 else 12; t0 = time.time()
ix = ct.lnaKeyIndex(N); out = []
for key, mem in ix.items():
    if len(mem) < 2: continue
    cl = collections.defaultdict(list)
    for rl in mem:
        C = np.array(ct.lnaCartanMatrix(N, rl), dtype=int); S = C + C.T
        snf = smith_normal_form(sp.Matrix(S), domain=sp.ZZ); snf = tuple(sorted(abs(int(snf[i, i])) for i in range(N)))
        ev = np.linalg.eigvalsh(S.astype(float)); cl[(snf, (int((ev > 1e-9).sum()), int((ev < -1e-9).sum())))].append(ct.className(rl))
    if len(cl) > 1: out.append((len(mem), {k: sorted(v) for k, v in cl.items()}))
print('n=%d: %d keys, %d groups of >= 2 LNAs, %d split by (SNF, signature) [%.0fs]' % (N, len(ix), sum(len(m) > 1 for m in ix.values()), len(out), time.time() - t0))
for size, cl in sorted(out, key=lambda x: x[0]):
    print('  group of %d: %s' % (size, '; '.join('%d LNAs SNF tail %s sig %s root %s' % (len(v), k[0][-3:], k[1], v[0]) for k, v in cl.items())))
