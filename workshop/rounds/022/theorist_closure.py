"""T6 (round 022, theorist): for LNAs whose relations all have 2 arrows (digits in {0,2}), is every algebra reached within depth L (both
directions) a tree quiver (n-1 arrows) with only monomial relations of length 2? Prints the census of (arrows, relation shapes).
usage: theorist_closure.py N L   (repository root)."""
import sys, os, copy, itertools, collections
sys.path.insert(0, os.getcwd())
from quivermutation import coxeterTables as ct, nakayama as nk, pathAlgebra
from quivermutation import search as se, procedure as _pr
_pr.os = os
n, L = int(sys.argv[1]), int(sys.argv[2])
cen = collections.Counter(); seen = set()
for rl in sorted(ct.lnaStatus(n)):
    if max(rl) > 2: continue
    alg = nk.LinearNakayamaAlgebra(n, list(rl))
    def mk(dual):
        def v(a, p):
            if len(p) > L: return
            b = pathAlgebra.dualPathAlgebra(a) if dual else a
            ar = sorted(b.quiver.edges())
            shape = (len(ar), tuple(sorted((len(r), tuple(len(q) for q in r[:1])) for r in b.rels)))
            key = (tuple(ar), tuple(sorted(tuple(tuple(q) for q in r) for r in b.rels)))
            if key in seen: return
            seen.add(key)
            cen[(len(ar), tuple(sorted(set((len(r), tuple(sorted(len(q) for q in r))) for r in b.rels))))] += 1
        return v
    for dual, st in enumerate([alg, pathAlgebra.dualPathAlgebra(alg)]):
        se.mutationSearchDepthFirst(copy.deepcopy(st), L, [], 'lines', printOutput=False, visitor=mk(dual))
for k, c in sorted(cen.items()): print(k, c)
print("distinct labelled algebras", len(seen))
