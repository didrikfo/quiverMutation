"""Round 115 (theorist): HH^0,1,2 of every LNA and dual LNA of length n, grouped by Coxeter key (the key-guard class).  usage: theorist_lna.py N"""
import sys, time
sys.path.insert(0, '.'); N = int(sys.argv[1]); sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
exec(open('workshop/rounds/115/theorist_hh.py').read().split("if __name__")[0])
from collections import defaultdict
def hhof(alg):
    return hh(sorted(alg.quiver.nodes), sorted(alg.quiver.edges(keys=True)), procedure.relationsFrom(alg), top=2)
cls = defaultdict(list); t0 = time.time()
for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
    for tag, alg in (('F', lna), ('D', pathAlgebra.dualPathAlgebra(lna))):
        d = hhof(alg); cls[search._coxeterKeyOrNone(alg)].append((tag, d))
for k, L in sorted(cls.items(), key=lambda kv: len(kv[1])):
    vals = sorted({tuple(d[1]) if d[1] else None for _, d in L}, key=str)
    print('class size', len(L), 'dimA', sorted({d[0] for _, d in L}), 'distinct HH(0,1,2):', vals)
print('time %.0fs' % (time.time() - t0))
