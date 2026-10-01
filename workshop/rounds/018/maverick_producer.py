"""T6 (round 018, maverick): log the producing relation of every cord member (arrows >= n, no parallel arrows) found by the
walk from LNA index I at path length <= L (same visitor as workshop/rounds/015/toolsmith_cords.py, both directions).
For each member: relation kinds (monomial / sum), whether every arrow lying on an undirected cycle lies on a path of some
sum relation ("cycle covered by a sum relation"), and the shape (lengths of the two paths) of the sum relations.
usage: maverick_producer.py N L I [I ...]   (repository root)"""
import sys, os, collections, copy, time
sys.path.insert(0, os.getcwd())
import networkx as nx
from quivermutation import coxeterTables as ct, nakayama as nk, pathAlgebra
from quivermutation import search as se
n, L = int(sys.argv[1]), int(sys.argv[2]); idx = [int(x) for x in sys.argv[3:]]
def snap(a):
    ar = sorted(a.quiver.edges())
    if len(set(ar)) != len(ar): return None
    return (tuple(ar), tuple(sorted(tuple(tuple(p) for p in r) for r in a.rels)))
def members(alg):
    best = {}
    def mk(dual):
        def v(a, p):
            if len(p) > L: return
            s = snap(pathAlgebra.dualPathAlgebra(a) if dual else a)
            if s and len(s[0]) >= n and (s not in best or len(p) < best[s]): best[s] = len(p)
        return v
    for dual, st in enumerate([alg, pathAlgebra.dualPathAlgebra(alg)]):
        se.mutationSearchDepthFirst(copy.deepcopy(st), L, [], 'lines', printOutput=False, visitor=mk(dual))
    return best
def analyse(s):
    arrows, rels = s
    g = nx.Graph(); g.add_edges_from(arrows)
    cyc = set()
    for comp in nx.biconnected_components(g):
        h = g.subgraph(comp)
        if h.number_of_edges() >= 3:
            cyc |= {tuple(sorted(e)) for e in h.edges()}
    cyc2 = {tuple(sorted(e)) for e in cyc}
    cov = set(); kinds = []
    for r in rels:
        if len(r) == 1: kinds.append("mono%d" % (len(r[0]) - 1))
        else:
            kinds.append("sum" + "+".join(str(len(p) - 1) for p in r))
            for p in r:
                for a, b in zip(p, p[1:]): cov.add(tuple(sorted((a, b))))
    return kinds, cyc2 <= cov, len(arrows) - n + 1, cyc2 != set()
status = ct.lnaStatus(n); lst = sorted(status)
tot = collections.Counter(); shapes = collections.Counter()
for k in idx:
    rl = lst[k]; t = time.time()
    best = members(nk.LinearNakayamaAlgebra(n, list(rl)))
    c = collections.Counter(); sh = collections.Counter(); monoonly = 0
    for s, l in best.items():
        kinds, covered, cyclomatic, hascyc = analyse(s)
        c[("cycles", cyclomatic)] += 1
        c[("has sum rel", any(x.startswith("sum") for x in kinds))] += 1
        c[("all cycle arrows in a sum-relation support", covered)] += 1
        if not any(x.startswith("sum") for x in kinds): monoonly += 1
        for x in set(kinds):
            if x.startswith("sum"): sh[x] += 1
    print("lna", k, "".join(map(str, rl)), "members(L<=%d):" % L, len(best), dict(c), "monomial-only:", monoonly, "(%.0fs)" % (time.time() - t), flush=True)
    print("   sum-relation shapes (members containing one):", dict(sorted(sh.items())), flush=True)
    tot.update(c); shapes.update(sh)
print("TOTAL", dict(tot)); print("shapes", dict(sorted(shapes.items())))
