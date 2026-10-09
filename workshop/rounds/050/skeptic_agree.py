"""Round 050 (skeptic): does the independent test (Hom(T,T[-1]) = 0 = Hom(T,T[1]), T = T_v) agree with tiltingPlus + J? On all LNAs of length n, all duals,
and every vertex v with out-arrows: cross-tab (tiltingPlus, mine, J, cartan-congruence H = C(child)^T). Usage: skeptic_agree.py n"""
import sys
n = int(sys.argv[1]); sys.argv = ['x']
sys.path.insert(0, '.')
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import numpy as np
from collections import Counter
cnt = Counter()
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        rels = procedure.relationsFrom(alg)
        for v in sorted(alg.quiver.nodes):
            if not ap.arrowsOutOf(alg.quiver, v): continue
            tp = bool(tiltingPlus(alg.quiver, rels, v)); J = bool(perI(alg, v))
            V, res = stepTest(alg, v, -1)
            mine = not np.array(res[-1]).any() and not np.array(res[1]).any()
            gate = bool(mutation.mutationIsPossibleAtVertex(alg, v))
            cnt[(gate, tp, J, mine)] += 1
for k, c in sorted(cnt.items()): print('gate %s tiltingPlus %s J!=0 %s indep-tilt %s : %d' % (k[0], k[1], k[2], k[3], c))
