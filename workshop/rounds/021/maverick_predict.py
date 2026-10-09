"""T6 (round 021, maverick): which n LNAs have cord members (arrows >= n, no parallel arrows) within mutation-path length L?
Walk (both directions, same visitor as rounds/018/maverick_producer.py) from LNAs with index in [LO, HI], print per LNA
the relation lengths, the count of cord members at length <= L and the walk time, plus structural features.
usage: maverick_predict.py N L LO HI   (repository root); writes one line per LNA "idx seq ncord secs"."""
import sys, os, copy, time
sys.path.insert(0, os.getcwd())
from quivermutation import coxeterTables as ct, nakayama as nk, pathAlgebra
from quivermutation import search as se
n, L, lo, hi = [int(x) for x in sys.argv[1:5]]
lst = sorted(ct.lnaStatus(n))
def members(alg):
    best = {}
    def mk(dual):
        def v(a, p):
            if len(p) > L: return
            b = pathAlgebra.dualPathAlgebra(a) if dual else a
            if os.environ.get('MONO') and any(len(r) != 1 for r in b.rels): return  # MONO=1: monomial relations only
            ar = sorted(b.quiver.edges())
            if len(ar) < n or len(set(ar)) != len(ar): return
            s = (tuple(ar), tuple(sorted(tuple(tuple(q) for q in r) for r in b.rels)))
            if s not in best or len(p) < best[s]: best[s] = len(p)
        return v
    for dual, st in enumerate([alg, pathAlgebra.dualPathAlgebra(alg)]):
        se.mutationSearchDepthFirst(copy.deepcopy(st), L, [], 'lines', printOutput=False, visitor=mk(dual))
    return best
ks = [int(x) for x in os.environ['IDX'].split(',')] if os.environ.get('IDX') else range(lo, min(hi, len(lst) - 1) + 1)  # IDX=i,j,.. overrides LO HI
for k in ks:
    t = time.time(); b = members(nk.LinearNakayamaAlgebra(n, list(lst[k])))
    print(k, "".join(map(str, lst[k])), len(b), min(b.values()) if b else -1, "%.1f" % (time.time() - t), flush=True)
