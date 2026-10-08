"""Round 042 (theorist): explore J != 0 records of a round-041-style dump: C', G = C'^-1, Q coefficients.
Usage: theorist_explore.py dump.pkl"""
import sys, pickle
import numpy as np
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']
def prep(r):
    V, v, CA, CB = r['V'], r['v'], r['CA'], r['CB']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    return k, rm, rm @ CA @ rm.T
shown = 0; cnt = Counter()
for r in D['recs']:
    Jn = {i: j for i, j in r['J'].items() if j}
    if not Jn: continue
    k, rm, Cp = prep(r); i = r['V'].index(list(Jn)[0])
    G = np.rint(np.linalg.inv(Cp)).astype(int)
    cnt[(tuple(Cp[i]), )] += 0
    if shown < 3:
        shown += 1
        print('v', k, 'i', i, 'out', [r['V'].index(w) for w in r['out']]); print('CA\n', r['CA']); print("C'\n", Cp); print('CB-C\'\n', r['CB'] - Cp); print('G\n', G)
