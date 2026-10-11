"""T5 round 014: is the A5 algebra (E-080), padded by chains, in the Coxeter class of any LNA?
The Coxeter polynomial is invariant along legal mutation, and the guard keeps it fixed, so a guarded walk
from an LNA stays on that LNA's key. If the key of a hand-built algebra is on no LNA's key, no guarded walk
from an LNA at that n reaches it. Usage: scholar_a5key.py [nmax=9]"""
import sys; nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
sys.path.insert(0, '.'); sys.argv = ['x']
src = open('workshop/rounds/011/scholar_square.py').read().split("print(\"kind pre")[0]
exec(compile(src, 'sq', 'exec'))
from quivermutation import search
for n in range(5, nmax + 1):
    keys = {}
    for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
            k = search._coxeterKeyOrNone(alg)
            keys.setdefault(k, 0); keys[k] += 1
    print('n =', n, 'distinct LNA keys:', len(keys), flush=True)
    for pre in range(0, n - 4):
        post = n - 5 - pre
        alg, d = build(pre, post, 'long')
        k = search._coxeterKeyOrNone(alg)
        print('  A5 pre=%d post=%d key=%s in LNA keys: %s' % (pre, post, k, k in keys), flush=True)
