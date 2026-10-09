"""T6 (round 022, theorist): test the cord criterion 'cord member (arrows >= n, no parallel arrows) within depth L iff some relation has >= 3 arrows'
at any n, and report per-LNA minimum cord depth (up to L), plus the first cord member found.
usage: theorist_cordcrit.py N L [LO HI] [-v]   (repository root). Same visitor as rounds/021/maverick_predict.py."""
import sys, os, copy
sys.path.insert(0, os.getcwd())
from quivermutation import coxeterTables as ct, nakayama as nk, pathAlgebra
from quivermutation import search as se, procedure as _pr
_pr.os = os  # procedure.py (uncommitted edit by another worker) lacks `import os`
args = [a for a in sys.argv[1:] if a != '-v']; verbose = '-v' in sys.argv
n, L = int(args[0]), int(args[1])
lst = sorted(ct.lnaStatus(n)); lo = int(args[2]) if len(args) > 2 else 0; hi = int(args[3]) if len(args) > 3 else len(lst) - 1
def firstCord(alg):
    best = {}
    def mk(dual):
        def v(a, p):
            if len(p) > L: return
            b = pathAlgebra.dualPathAlgebra(a) if dual else a
            ar = sorted(b.quiver.edges())
            if len(ar) < n or len(set(ar)) != len(ar): return
            if not best or len(p) < best['d']:
                best.update(d=len(p), p=list(p), ar=ar, rels=[[list(q) for q in r] for r in b.rels], dual=dual)
        return v
    for dual, st in enumerate([alg, pathAlgebra.dualPathAlgebra(alg)]):
        se.mutationSearchDepthFirst(copy.deepcopy(st), L, [], 'lines', printOutput=False, visitor=mk(dual))
    return best
mis = 0; cnt = 0
ks = [i for i, r in enumerate(lst) if "".join(map(str, r)) in os.environ['NAMES'].split(',')] if os.environ.get('NAMES') else range(lo, min(hi, len(lst) - 1) + 1)  # NAMES=d1,d2,.. restricts to those LNAs
for k in ks:
    rl = lst[k]; b = firstCord(nk.LinearNakayamaAlgebra(n, list(rl)))
    pred = max(rl) >= 3; got = bool(b); cnt += 1
    if pred != got: mis += 1
    print(k, "".join(map(str, rl)), "pred", int(pred), "found", int(got), "depth", b.get('d', -1), ("" if not (verbose and b) else (b['p'], b['dual'], b['ar'], b['rels'])), flush=True)
print("n", n, "L", L, "LNAs", cnt, "mismatches", mis)
