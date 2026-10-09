"""Round 042 (theorist): for the J != 0 steps with u = e_w - e_i, find which coincidences F^k e_a = +-F^l e_b (a, b in {i, w}; |k|,|l| <= 3;
F = Z Z^-T, Z = C_A restricted to V\\v) hold, and on how many records (split into 'type 1': F e_w = e_i, and the rest).
Usage: theorist_orbits.py dump.pkl"""
import sys, pickle, itertools
import numpy as np
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']; m = n - 1
tot = Counter(); recs = []
for r in D['recs']:
    Jn = {i: j for i, j in r['J'].items() if j}
    if not Jn or len(r['out']) != 1: continue
    V, v = r['V'], r['v']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    Cp = rm @ r['CA'] @ rm.T; keep = [a for a in range(n) if a != k]
    Z = Cp[np.ix_(keep, keep)]; w = keep.index(V.index(r['out'][0])); i = keep.index(V.index(list(Jn)[0]))
    Y = np.rint(np.linalg.inv(Z)).astype(int); F = Z @ Y.T; Fi = np.rint(np.linalg.inv(F)).astype(int)
    recs.append((F, Fi, w, i))
def pw(F, Fi, k): 
    M = np.eye(m, dtype=int)
    for _ in range(abs(k)): M = M @ (F if k > 0 else Fi)
    return M
for typ in ('type1', 'other'):
    c = Counter(); N = 0
    for F, Fi, w, i in recs:
        t1 = bool((F[:, w] == np.eye(m, dtype=int)[i]).all())
        if (typ == 'type1') != t1: continue
        N += 1
        vecs = {(a, kk): pw(F, Fi, kk)[:, idx] for (a, idx) in (('i', i), ('w', w)) for kk in range(-3, 4)}
        keys = list(vecs)
        for p, q in itertools.combinations(keys, 2):
            if p[0] == q[0] and p[1] == q[1]: continue
            if (vecs[p] == vecs[q]).all(): c[(p, '=', q)] += 1
            elif (vecs[p] == -vecs[q]).all(): c[(p, '=-', q)] += 1
    print(typ, 'records', N)
    for key, cc in sorted(c.items(), key=lambda t: -t[1]):
        if cc == N: print('   always:', key)
