"""H-017 control at path length L >= 5 (extends rounds/007/maverick_control.py) plus shortest-path check.
For LNAs of order N (status from coxeterTables): walk to DEPTH (reachedQuipuAlgebras, exhaustive DFS,
no Visited pruning, so the recorded path is the shortest within DEPTH). Members with recorded path
length == LEN (hereditary excluded) are taken PER per LNA; for each:
  shortest: rewalk the source at depth LEN-1 and confirm the certificate is NOT reached (and at LEN that it is);
  control:  rebuild the member as a path algebra, verify-style search (linesReachedFrom) at depth LEN
            returns the source LNA (hit), and at depth LEN-1 does not (neg).
usage: maverick_control5.py N DEPTH LEN PER [LIMIT_LNAS] [--plan]"""
import sys, time, collections
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk, search as se, pathAlgebra, lines
args = [a for a in sys.argv[1:] if not a.startswith("--")]; plan = "--plan" in sys.argv
n, depth, L, per = map(int, args[:4]); lim = int(args[4]) if len(args) > 4 else 10**9
def algOf(cert):
    a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, n + 1))
    for t, h in cert[0]: a.add_arrow(t + 1, h + 1)
    for p in cert[1]: a.add_rel([[v + 1 for v in p]])
    return a
def names(alg, d):
    return {"".join(str(x) for x in lines.relationStringToLineRelLengths(n, s)) for s in se.linesReachedFrom(alg, d)}
status = ct.lnaStatus(n); T = collections.Counter(); T0 = time.time(); done = 0
for rl in sorted(status):
    if done >= lim: break
    src = "".join(map(str, rl)); alg = nk.LinearNakayamaAlgebra(n, list(rl)); t = time.time()
    reached = qr.reachedQuipuAlgebras(alg, depth)
    lens = collections.Counter(len(p) for c, p in reached.items() if not all(h == a + 1 for a, h in c[0]))
    mem = [(c, p) for c, p in reached.items() if len(p) == L and not all(h == a + 1 for a, h in c[0])]
    mem.sort(key=lambda x: -len(x[0][1]))
    if plan: print(src, "walk %.1fs" % (time.time() - t), "members by path length", dict(sorted(lens.items())), flush=True)
    if not mem: T["no member at L"] += 1; done += 1; continue
    done += 1
    if plan: continue
    below = qr.reachedQuipuAlgebras(alg, L - 1)
    for c, p in mem[:per]:
        T["short_ok" if c not in below else "SHORT_FAIL"] += 1
        a = algOf(c)
        hit = src in names(a, L); neg = src not in names(a, L - 1)
        T["hit" if hit else "MISS"] += 1; T["neg_ok" if neg else "NEG_FAIL"] += 1
        print(src, "L", L, "rels", len(c[1]), "hit", hit, "neg", neg, "%.1fs" % (time.time() - t), flush=True)
print(dict(T), "total %.0fs" % (time.time() - T0))
