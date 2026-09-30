"""H-017 walk: (cords, relations) of quipu-with-relations algebras PROVED in the class
of each LNA outside a quipu class, by mutation walk to depth D.
usage: python workshop/rounds/004/maverick_reached.py N DEPTH [CLASSNAME ...]   (default: every LNA outside a quipu class)"""
import sys, collections, networkx as nx
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk, quipuForms as qf
n, depth = int(sys.argv[1]), int(sys.argv[2])
status = ct.lnaStatus(n)
outside = sorted(r for r, v in status.items() if v != ct.QUIPU and (len(sys.argv) < 4 or ct.className(r) in sys.argv[3:]))
allpairs = collections.Counter(); worst = []
for rl in outside:
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
