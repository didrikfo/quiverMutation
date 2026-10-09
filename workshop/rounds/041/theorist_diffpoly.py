"""Round 041 (theorist): the polynomial Q(x) = det(xC_B + C_B^T) - det(xC' + C'^T) for J != 0 records of a dump. Prints coefficient vectors (x^0..x^n) with counts.
Usage: theorist_diffpoly.py dump.pkl"""
import sys, pickle
import numpy as np
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']
xs = np.arange(0, n + 1); Vd = np.vander(xs.astype(float), n + 1, increasing=True)
def coef(C): return np.linalg.solve(Vd, np.array([np.linalg.det(t * C + C.T) for t in xs]))
c = Counter(); lowest = Counter()
for r in D['recs']:
    Jn = {i: j for i, j in r['J'].items() if j}
    if not Jn: continue
    V, v, CA, CB = r['V'], r['v'], r['CA'], r['CB']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    Cp = rm @ CA @ rm.T
    q = tuple(int(round(a)) for a in coef(CB) - coef(Cp)); c[(q, tuple(sorted(Jn.values())))] += 1
    lowest[next(i for i, a in enumerate(q) if a)] += 1
print('n', n, 'J != 0 records', sum(c.values()), ' lowest nonzero degree of Q:', dict(lowest))
for (q, j), m in sorted(c.items(), key=lambda t: -t[1])[:12]: print('  Q =', q, ' dimJ', j, ' count', m)
