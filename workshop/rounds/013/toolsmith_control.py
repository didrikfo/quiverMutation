"""Depth-6 positive control for the H-017 mutation search, with node counts (round 013, toolsmith).
For LNAs at order N: walk out to WALKDEPTH (qr.reachedQuipuAlgebras), take hereditary-or-relations members whose recorded
path has length exactly L, rebuild each as a path algebra and search it at depth SEARCH (default L) -> must return the source LNA;
with SHORT=1 search at L-1 (negative control).  Prints per member: path length, nodes, distinct states, found, seconds.
usage: toolsmith_control.py N WALKDEPTH L [PER_LNA] [FIRST_LNA_INDEX] [LAST_LNA_INDEX]   (run from the repository root)"""
import sys, os, time, collections, importlib.util
sys.path.insert(0, os.getcwd())
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk, pathAlgebra
import copy
from quivermutation import search as se, lines, fingerprint
n, wd, L = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
per = int(sys.argv[4]) if len(sys.argv) > 4 else 1
lo = int(sys.argv[5]) if len(sys.argv) > 5 else 0
hi = int(sys.argv[6]) if len(sys.argv) > 6 else 10**9
sd = L - 1 if os.environ.get("SHORT") else L
def algOf(cert):
    a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, n + 1))
    for t, h in cert[0]: a.add_arrow(t + 1, h + 1)
    for p in cert[1]: a.add_rel([[v + 1 for v in p]])
    return a
def counted(alg, d):
    visits = [0]; keys = set()
    def visit(a, p): visits[0] += 1; keys.add(fingerprint.canonicalKey(a))
    names = set()
    for i, st in enumerate([alg, pathAlgebra.dualPathAlgebra(alg)]):
        col = []
        se.mutationSearchDepthFirst(copy.deepcopy(st), d, [], 'lines', printOutput=False, collected=col, visitor=visit)
        for f in lines.mutationListLineCleanup(col, printOutput=False):
            a = f[0] if i == 0 else pathAlgebra.dualPathAlgebra(f[0])
            names.add("".join(str(x) for x in lines.relationStringToLineRelLengths(n, lines.relSetToString(sorted(se._renumberedLine(a).rels)))))
    return names, visits[0], len(keys)
status = ct.lnaStatus(n); stat = collections.Counter(); tot = [0, 0]
for k, rl in enumerate(sorted(status)):
    if k < lo or k > hi: continue
    src = "".join(str(x) for x in rl)
    t = time.time(); reached = qr.reachedQuipuAlgebras(nk.LinearNakayamaAlgebra(n, list(rl)), wd)
    mem = [(len(p), c) for c, p in reached.items() if len(p) == L and not all(h == t_ + 1 for t_, h in c[0])]
    mem.sort(key=lambda x: len(x[1][1]))
    print("lna", k, src, "members with L =", L, ":", len(mem), "(walk %.0fs)" % (time.time() - t), flush=True)
    for l, c in mem[:per]:
        t = time.time(); got, nv, nd = counted(algOf(c), sd)
        stat[src in got] += 1; tot[0] += nv; tot[1] += nd
        print("  lna", k, src, "L", l, "search", sd, "rels", len(c[1]), "found", src in got, "nodes", nv, "distinct", nd, "%.0fs" % (time.time() - t), flush=True)
print("found counts:", dict(stat), "total nodes", tot[0], "distinct (summed)", tot[1])
