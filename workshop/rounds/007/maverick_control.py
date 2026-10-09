"""Positive control for the H-017 mutation search (maverick_verify.py / families.verify).
For LNAs at order N: walk out to DEPTH (qr.reachedQuipuAlgebras), take quipu-with-relations members
with a recorded path of length L, then run the verify-style search (se.linesReachedFrom, both directions)
FROM the member at depth D and ask whether the source LNA (and any LNA of its class) comes back.
SHORT=1 in the environment searches one level shallower than the path (a negative control).
usage: maverick_control.py N DEPTH SEARCHDEPTH [PER_LNA] [CLASSNAME ...]"""
import sys, os, collections, time
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk, search as se, pathAlgebra, lines
n, depth, sdepth = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
per = int(sys.argv[4]) if len(sys.argv) > 4 else 3
only = sys.argv[5:]
def algOf(cert):
    a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, n + 1))
    for t, h in cert[0]: a.add_arrow(t + 1, h + 1)
    for p in cert[1]: a.add_rel([[v + 1 for v in p]])
    return a
def names(alg, d):
    red = se.linesReachedFrom(alg, d)
    return {"".join(str(x) for x in lines.relationStringToLineRelLengths(n, s)) for s in red}
status = ct.lnaStatus(n)
stat = collections.Counter(); bylen = collections.defaultdict(collections.Counter)
for rl in sorted(status):
    nm = ct.className(rl)
    if only and nm not in only: continue
    src = "".join(str(x) for x in rl)
    alg = nk.LinearNakayamaAlgebra(n, list(rl))
    reached = qr.reachedQuipuAlgebras(alg, depth)
    mem = [(len(p), c, p) for c, p in reached.items() if not all(h == t + 1 for t, h in c[0])]
    mem.sort(key=lambda x: (-x[0], len(x[1][1])))
    for L, c, p in mem[:per]:
        t = time.time(); got = names(algOf(c), (L - 1) if os.environ.get("SHORT") else max(L, sdepth))
        hit = src in got; cls = bool(got & {"".join(str(x) for x in r) for r in status if ct.className(r) == nm})
        stat[(hit, cls)] += 1; bylen[L][hit] += 1
        print(nm, "L", L, "rels", len(c[1]), "source found", hit, "| LNAs reached", len(got), "%.1fs" % (time.time() - t), flush=True)
print("(source found, class found):", dict(stat)); print("by path length:", {k: dict(v) for k, v in sorted(bylen.items())})
