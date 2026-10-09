"""T6 (round 018, maverick): code-level positive control for the MONO=1 filter of workshop/rounds/015/toolsmith_cords.py.
Seeds the same visitor with a hand-built MONOMIAL cord algebra (A3 triangle 1->2->3, 1->3, zero relation 1->2->3, plus leaf
3->4 at n = 4) instead of an LNA, so the visitor must return it at depth 0 under MONO=1; then lists what its orbit holds
(monomial cord members by path length), to see whether the filter ever returns a *different* one.
usage: MONO=1 maverick_filtercheck.py L   (repository root)"""
import sys, os, copy, collections
sys.path.insert(0, os.getcwd())
from quivermutation import pathAlgebra, search as se
L = int(sys.argv[1]); n = 4
a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, n + 1))
for t, h in [(1, 2), (2, 3), (1, 3), (3, 4)]: a.add_arrow(t, h)
a.add_rel([[1, 2, 3]])
def snap(x):
    if os.environ.get("MONO") and any(len(r) != 1 for r in x.rels): return None
    ar = sorted(x.quiver.edges())
    if len(set(ar)) != len(ar): return None
    return (tuple(ar), tuple(sorted(tuple(tuple(p) for p in r) for r in x.rels)))
best = {}
def v(x, p):
    if len(p) > L: return
    s = snap(x)
    if s and len(s[0]) >= n and (s not in best or len(p) < best[s]): best[s] = len(p)
for st in [a, pathAlgebra.dualPathAlgebra(a)]:
    se.mutationSearchDepthFirst(copy.deepcopy(st), L, [], 'lines', printOutput=False, visitor=v)
print("MONO", bool(os.environ.get("MONO")), "L", L, "cord members", len(best), "by path length", dict(sorted(collections.Counter(best.values()).items())))
for s, l in sorted(best.items(), key=lambda t: t[1])[:8]: print(" ", l, s)
