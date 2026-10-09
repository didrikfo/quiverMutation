"""Round 046 (skeptic), response to referee item 1: control.  For key classes of LNAs (and duals) at length n, search P in GL_n(Z) (box B) with
P C_X P^T = C_Y for ALL pairs of distinct Cartan matrices in the class.  Equal key => P expected.  Usage: skeptic_pairs.py n B class_indices..."""
import sys, time, itertools
ARGV = sys.argv; sys.argv = ['x']
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/046')
import numpy as np
from quivermutation import nakayama as nk, pathAlgebra, search, invariants as qi
src = open('workshop/rounds/046/skeptic_iso.py').read()
exec(src[src.index('def search('):src.index('for k, r in')].replace('def search(', 'def psearch('))
n = int(ARGV[1]); B = int(ARGV[2]); want = [int(x) for x in ARGV[3:]]
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
for ci in want:
    mats = {}
    for alg in classes[order[ci]]:
        C = [[int(x) for x in r] for r in qi.cartanMatrix(alg).tolist()]
        mats[str(C)] = C
    ms = list(mats.values()); found = none = to = 0; t0 = time.time()
    for X, Y in itertools.combinations(ms, 2):
        res, _ = psearch(X, Y, B)
        if res is None: to += 1
        elif res is False: none += 1
        else: found += 1
    print('class', ci, 'algebras', len(classes[order[ci]]), 'distinct Cartan', len(ms), 'pairs', found + none + to, 'P found', found, 'none', none, 'timeout', to, '%ds' % (time.time() - t0), flush=True)
