"""Round 046 (skeptic): follow-up to skeptic_iso.py: the one n = 7 c1 failing step (index 12) with no P in box 1, box 2 and 3; and a negative control
(random pair, equal charpoly, separated by skeptic_inv: the search must find nothing).  Usage: skeptic_iso12.py c1.pkl"""
import sys, pickle, random
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/046'); ARGV = sys.argv; sys.argv = ['x', ARGV[1], '1']
import numpy as np
src = open('workshop/rounds/046/skeptic_iso.py').read().split("for k, r in enumerate")[0]; exec(compile(src, 'iso', 'exec'))
r = [x for x in recs if x['kind'] == 'fail'][12]
for b in (1, 2, 3):
    for name, CB in (('C_B', r['CB']), ('C_B^T', np.array(r['CB']).T.tolist())):
        res, nc = search(r['CA'], CB, b, 300); print('box', b, name, 'timeout' if res is None else ('FOUND' if res is not False else 'none'), nc, flush=True)
        if res is not None and res is not False: print(res.tolist()); break
import skeptic_inv as si
random.seed(1); n = 7
def rnd():
    C = np.eye(n, dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.35: C[i, j] = random.choice([1, 1, 1, 2, -1])
    return C
groups = {}; done = False
while not done:
    C = rnd(); inv = si.invariants(C.tolist(), mmax=3); k = inv['charpoly']
    for (C0, i0) in groups.get(k, []):
        if si.diff(inv, i0):
            res, nc = search(C0.tolist(), C.tolist(), 1, 120); print('neg control: invariants differ in', si.diff(inv, i0)[:2], '-> iso search', 'timeout' if res is None else ('FOUND' if res is not False else 'none')); done = True; break
    groups.setdefault(k, []).append((C, inv))
