"""Round 042 (theorist): check C_B = r C_A r^T + e_v J^T (E-136) on every gate-admitted record of a dump (closes the round-041 review item for n = 7).
Usage: theorist_cbcheck.py dump.pkl"""
import sys, pickle
import numpy as np
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']; c = Counter()
for r in D['recs']:
    V, v = r['V'], r['v']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    H = np.zeros((n, n), dtype=int)
    for i, j in r['J'].items(): H[k, V.index(i)] = j
    ok = bool((r['CB'] == rm @ r['CA'] @ rm.T + H).all()); c[ok, bool(any(r['J'].values()))] += 1
print('n', n, '(C_B = C\'+H holds?, J != 0):', dict(c))
