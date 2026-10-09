"""H-017 control at n = 8 with a NON-HEREDITARY start (round 014, maverick). Extends 013/toolsmith_control.py.
For LNA k at order N: walk to WALKDEPTH (qr.reachedQuipuAlgebras), keep members whose recorded path has length exactly L
AND which have >= MINREL relations (non-hereditary), take up to PER with fewest relations (ties: fewest arrows' cords),
rebuild, search at depth L (SHORT=1: L-1) with node counts. --plan: list members only, no search.
usage: maverick_control8.py N WALKDEPTH L MINREL PER FIRST LAST [--plan]   (repository root)"""
import sys, os, time, collections, copy
sys.path.insert(0, os.getcwd())
plan = "--plan" in sys.argv; sys.argv = [a for a in sys.argv if a != "--plan"]
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk, pathAlgebra
from quivermutation import search as se, lines, fingerprint
n, wd, L, minrel, per, lo, hi = [int(x) for x in sys.argv[1:8]]
sd = L - 1 if os.environ.get("SHORT") else L
def algOf(cert):
    a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, n + 1))
    for t, h in cert[0]: a.add_arrow(t + 1, h + 1)
    for p in cert[1]: a.add_rel([[v + 1 for v in p]])
    return a
def counted(alg, d):
    visits = [0]; keys = set(); names = set()
    def visit(a, p): visits[0] += 1; keys.add(fingerprint.canonicalKey(a))
    for i, st in enumerate([alg, pathAlgebra.dualPathAlgebra(alg)]):
        col = []
        se.mutationSearchDepthFirst(copy.deepcopy(st), d, [], 'lines', printOutput=False, collected=col, visitor=visit)
        for f in lines.mutationListLineCleanup(col, printOutput=False):
            a = f[0] if i == 0 else pathAlgebra.dualPathAlgebra(f[0])
            names.add("".join(str(x) for x in lines.relationStringToLineRelLengths(n, lines.relSetToString(sorted(se._renumberedLine(a).rels)))))
    return names, visits[0], len(keys)
status = ct.lnaStatus(n); print("n", n, "LNAs", len(status), flush=True)
stat = collections.Counter(); tot = [0, 0]
for k, rl in enumerate(sorted(status)):
    if k < lo or k > hi: continue
    src = "".join(str(x) for x in rl)
    t = time.time(); reached = qr.reachedQuipuAlgebras(nk.LinearNakayamaAlgebra(n, list(rl)), wd)
    mem = [(len(p), c) for c, p in reached.items() if len(p) == L and len(c[1]) >= minrel]
    mem.sort(key=lambda x: ((-1 if os.environ.get("HIGH") else 1) * len(x[1][1]), len(x[1][0])))   # HIGH=1: most relations first
    hist = collections.Counter(len(c[1]) for _, c in mem)
    print("lna", k, src, "members L =", L, "rels>=%d:" % minrel, len(mem), "rels hist", dict(sorted(hist.items())), "(walk %.0fs)" % (time.time() - t), flush=True)
    if plan: continue
    for l, c in mem[:per]:
        t = time.time(); got, nv, nd = counted(algOf(c), sd)
        stat[src in got] += 1; tot[0] += nv; tot[1] += nd
        print("  lna", k, src, "L", l, "search", sd, "rels", len(c[1]), "arrows", len(c[0]), "found", src in got, "nodes", nv, "distinct", nd, "%.0fs" % (time.time() - t), flush=True)
print("found counts:", dict(stat), "total nodes", tot[0], "distinct (summed)", tot[1])
