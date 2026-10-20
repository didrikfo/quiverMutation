"""H-017 defect-1 source control (E-185 left open). Defect = relations - cords.
Forward-walk each named LNA outside a quipu class at n=9 to DEPTH, list members of path length L with defect <= D.
Without --run: plan only (counts, times). With --run: rebuild up to PER such members as path algebras,
search back at depth L and L-1, report hit / negative.
usage: LO=start maverick_defect1_rest.py N DEPTH L D PER NAME [--run]; run members sel[LO:LO+PER] (sharded copy of round 069 maverick_defect1.py)"""
import sys, time, collections, networkx as nx
from quivermutation import quipuRelations as qr, nakayama as nk, search as se, pathAlgebra, lines, quipuForms as qf
args = [a for a in sys.argv[1:] if not a.startswith("--")]; run = "--run" in sys.argv
n, depth, L, D, per = map(int, args[:5]); names_ = args[5:]
import os; LO = int(os.environ.get("LO", "0"))
def algOf(cert):
    a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, n + 1))
    for t, h in cert[0]: a.add_arrow(t + 1, h + 1)
    for p in cert[1]: a.add_rel([[v + 1 for v in p]])
    return a
def reach(alg, d):
    return {"".join(str(x) for x in lines.relationStringToLineRelLengths(n, s)) for s in se.linesReachedFrom(alg, d)}
def cr(c):
    g = nx.Graph(); g.add_edges_from(c[0]); par = qf.quipuParameters(g)
    return sum(1 for x in par[1] if x > 0), len(c[1])
T = collections.Counter(); T0 = time.time()
for src in names_:
    t = time.time(); alg = nk.LinearNakayamaAlgebra(n, [int(x) for x in src])
    reached = qr.reachedQuipuAlgebras(alg, depth)
    mem = [(c, p, cr(c)) for c, p in reached.items() if not all(h == a + 1 for a, h in c[0])]
    sel = [(c, p, kr) for c, p, kr in mem if len(p) == L and kr[1] - kr[0] <= D and kr[0] >= 1]
    print(src, "members", len(mem), "L=%d cords>=1 defect<=%d: %d" % (L, D, len(sel)),
          collections.Counter(kr for _, _, kr in sel), "%.0fs" % (time.time() - t), flush=True)
    for c, p, kr in sel:
        print("   member", src, "path", p, "arrows", sorted(c[0]), "rels", [list(x) for x in c[1]], "cells", kr, flush=True)
    if not run: continue
    for c, p, (k, r) in sel[LO:LO + per]:
        a = algOf(c); hit = src in reach(a, L); neg = src not in reach(a, L - 1)
        T["hit" if hit else "MISS"] += 1; T["neg_ok" if neg else "NEG_FAIL"] += 1
        print("  ", src, "cords", k, "rels", r, "hit", hit, "neg", neg, "%.0fs" % (time.time() - t), flush=True)
print(dict(T), "total %.0fs" % (time.time() - T0))
