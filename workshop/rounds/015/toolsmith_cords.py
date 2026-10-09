"""H-017 control at n = 8 whose start has CORDS (arrows >= n) and relations >= MINREL (round 015, toolsmith).
reachedQuipuAlgebras cannot supply one (certificate() keeps quipu trees; every member at n = 8, LNA 0 has 7 arrows),
so members are collected here by a visitor on the raw search (both directions): any visited quiver with no parallel arrows,
relations monomial or sums (MONO=1: monomial only; none exist at n=7,8 within depth 5-6, see 015 toolsmith.md), arrows >= n, rels >= MINREL, recorded at mutation-path length exactly L (smallest seen for that labelled
algebra). Pick PER of them (most arrows, then fewest relations; HIGH=1 most relations; ties by key), rebuild, search from the
member at depth L (and L-1 with SHORT=1 or always with BOTH=1), with node counts. --plan: histogram only.
usage: toolsmith_cords.py N L MINREL PER FIRST LAST [--plan]   (repository root)"""
import sys, os, time, collections, copy
sys.path.insert(0, os.getcwd())
plan = "--plan" in sys.argv; sys.argv = [a for a in sys.argv if a != "--plan"]
from quivermutation import coxeterTables as ct, nakayama as nk, pathAlgebra
from quivermutation import search as se, lines, fingerprint
n, L, minrel, per, lo, hi = [int(x) for x in sys.argv[1:7]]
def algOf(arrows, rels):
    a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, n + 1))
    for t, h in arrows: a.add_arrow(t, h)
    for r in rels: a.add_rel([list(p) for p in r])
    return a
def snap(a):
    if os.environ.get("MONO") and any(len(r) != 1 for r in a.rels): return None
    ar = sorted(a.quiver.edges())
    if len(set(ar)) != len(ar): return None
    return (tuple(ar), tuple(sorted(tuple(tuple(p) for p in r) for r in a.rels)))
def members(alg):
    best = {}
    def mk(dual):
        def v(a, p):
            if len(p) > L: return
            s = snap(pathAlgebra.dualPathAlgebra(a) if dual else a)
            if s and len(s[0]) >= n and len(s[1]) >= minrel and (s not in best or len(p) < best[s]): best[s] = len(p)
        return v
    for dual, st in enumerate([alg, pathAlgebra.dualPathAlgebra(alg)]):
        se.mutationSearchDepthFirst(copy.deepcopy(st), L, [], 'lines', printOutput=False, visitor=mk(dual))
    return [s for s, l in best.items() if l == L], best
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
    t = time.time(); mem, best = members(nk.LinearNakayamaAlgebra(n, list(rl)))
    hist = collections.Counter((len(s[0]), len(s[1])) for s in mem)
    allh = collections.Counter((l, len(s[0]), len(s[1])) for s, l in best.items())
    print("lna", k, src, "cord members at L =", L, "rels>=%d:" % minrel, len(mem), "(arrows,rels) hist", dict(sorted(hist.items())), "(walk %.0fs)" % (time.time() - t), flush=True)
    print("  all cord members by (L,arrows,rels):", dict(sorted(allh.items())), flush=True)
    if plan: continue
    mem.sort(key=lambda s: (-len(s[0]), (-1 if os.environ.get("HIGH") else 1) * len(s[1]), s))
    for s in mem[:per]:
        for d in ([L, L - 1] if os.environ.get("BOTH") else [L - 1 if os.environ.get("SHORT") else L]):
            t = time.time(); got, nv, nd = counted(algOf(*s), d)
            stat[(d == L, src in got)] += 1; tot[0] += nv; tot[1] += nd
            print("  lna", k, src, "search", d, "arrows", len(s[0]), "rels", len(s[1]), "found", src in got, "nodes", nv, "distinct", nd, "%.0fs" % (time.time() - t), "member", s, flush=True)
print("(at-depth, found) counts:", dict(stat), "total nodes", tot[0], "distinct (summed)", tot[1])
