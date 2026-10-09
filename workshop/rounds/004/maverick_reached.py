"""H-017 walk: (cords, relations) of quipu-with-relations algebras PROVED in the class
of each LNA outside a quipu class, by mutation walk to depth D.
usage: python workshop/rounds/004/maverick_reached.py N DEPTH [CLASSNAME ...] [--budget-hours H]   (default: every LNA outside a quipu class)
With --budget-hours, stops between LNAs once H hours are spent and exits 2 (partial tallies printed)."""
import sys, time, collections, networkx as nx
_t0 = time.time(); _budget = None
if "--budget-hours" in sys.argv:
    _i = sys.argv.index("--budget-hours"); _budget = float(sys.argv[_i + 1]) * 3600; del sys.argv[_i:_i + 2]
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk, quipuForms as qf
n, depth = int(sys.argv[1]), int(sys.argv[2])
status = ct.lnaStatus(n)
outside = sorted(r for r, v in status.items() if v != ct.QUIPU and (len(sys.argv) < 4 or ct.className(r) in sys.argv[3:]))
allpairs = collections.Counter(); worst = []
_partial = False
for rl in outside:
    if _budget is not None and time.time() - _t0 > _budget:
        _partial = True; break
    alg = nk.LinearNakayamaAlgebra(n, list(rl))
    reached = qr.reachedQuipuAlgebras(alg, depth)
    pairs = collections.Counter()
    for cert, path in reached.items():
        if all(h == t + 1 for t, h in cert[0]): continue
        g = nx.Graph(); g.add_edges_from(cert[0])
        par = qf.quipuParameters(g)
        cords = sum(1 for x in par[1] if x > 0)
        r = len(cert[1]); pairs[(cords, r)] += 1; allpairs[(cords, r)] += 1
        if r <= cords: worst.append((ct.className(rl), qr.describeCertificate(cert), path))
    print(ct.className(rl), len(reached), "min(rel-cords) =", min((r-c for c, r in pairs), default=None),
          sorted(pairs), flush=True)
print("ALL pairs", sorted(allpairs.items()))
print("below/on diagonal:", len(worst))
for w in worst[:10]: print(w)
if _partial:
    print("BUDGET SPENT: partial"); sys.exit(2)
