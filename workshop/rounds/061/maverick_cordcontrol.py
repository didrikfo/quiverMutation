"""H-017 positive control restricted to tree quipus WITH cords (>=1 cord) and monomial relations.
For each LNA outside a quipu class at order N: forward walk to DEPTH (qr.reachedQuipuAlgebras), tally
(cords, relations) by recorded path length; with --run, take up to PER members of path length LEN with
cords >= 1, rebuild as a path algebra, search back (search.linesReachedFrom) at depth LEN and LEN-1.
usage: maverick_cordcontrol.py N DEPTH LEN PER [LIMIT_LNAS] [--run] [--onlyname NAME]"""
import sys, time, collections, networkx as nx
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk, search as se, pathAlgebra, lines, quipuForms as qf
args = [a for a in sys.argv[1:] if not a.startswith("--")]; run = "--run" in sys.argv
n, depth, L, per = map(int, args[:4]); lim = int(args[4]) if len(args) > 4 else 10**9
def algOf(cert):
    a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, n + 1))
    for t, h in cert[0]: a.add_arrow(t + 1, h + 1)
    for p in cert[1]: a.add_rel([[v + 1 for v in p]])
    return a
def names(alg, d):
    return {"".join(str(x) for x in lines.relationStringToLineRelLengths(n, s)) for s in se.linesReachedFrom(alg, d)}
def cr(c):
    g = nx.Graph(); g.add_edges_from(c[0]); par = qf.quipuParameters(g)
    return sum(1 for x in par[1] if x > 0), len(c[1])
status = ct.lnaStatus(n); outside = sorted(r for r, v in status.items() if v != ct.QUIPU)
T = collections.Counter(); T0 = time.time(); done = 0; table = collections.Counter()
for rl in outside:
    if done >= lim: break
    done += 1
    src = "".join(map(str, rl)); alg = nk.LinearNakayamaAlgebra(n, list(rl)); t = time.time()
    reached = qr.reachedQuipuAlgebras(alg, depth)
    mem = [(c, p, cr(c)) for c, p in reached.items() if not all(h == a + 1 for a, h in c[0])]
    for c, p, (k, r) in mem: table[(len(p), k, r)] += 1
    sel = [(c, p, kr) for c, p, kr in mem if len(p) == L and kr[0] >= 1]
    sel.sort(key=lambda x: (x[2][1] - x[2][0], -x[2][0]))   # nearest the diagonal first
    print(src, "members", len(mem), "L=%d cords>=1: %d" % (L, len(sel)), "%.1fs" % (time.time() - t), flush=True)
    if not run: continue
    for c, p, (k, r) in sel[:per]:
        a = algOf(c); hit = src in names(a, L); neg = src not in names(a, L - 1)
        T["hit" if hit else "MISS"] += 1; T["neg_ok" if neg else "NEG_FAIL"] += 1
        print("  ", src, "cords", k, "rels", r, "hit", hit, "neg", neg, "%.0fs" % (time.time() - t), flush=True)
print("table (pathlen,cords,rels):count"); 
for key in sorted(table): print(key, table[key])
print(dict(T), "total %.0fs" % (time.time() - T0))
