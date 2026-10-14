"""Item 4 of the round-065 review: Cartan check on the pairs with excess r7 > g (theorist).

Same sample as theorist_g2exact.py.  For every (quiver, i) with r7 > g and gate-admitted (J = 0, so
the step is tilting and R C R^T = Cartan(End T) must hold), compare R C R^T with the Cartan matrix
of (a) the raw rewrite (mutateAtVertex: all step 7 relations, including the excess rbar.s ones) and
(b) the rewrite after pr._minimalGenerators (excess dropped).  Also counts relations out of i*.
Child with an oriented cycle is reported separately (vertex-simple Cartan unreliable there).
Run: timeout 10m .venv/bin/python workshop/rounds/065/theorist_g2cartan.py [seeds] [per_seed]
"""
import sys, random
sys.path.insert(0, "."); sys.path.insert(0, "workshop/rounds/065")
import networkx as nx
import theorist_g2exact as G
from quivermutation import procedure as pr, arrowPaths as ap

seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 12
per = int(sys.argv[2]) if len(sys.argv) > 2 else 250
t = dict(excess=0, admitted=0, refused=0, cyc=0, rawOK=0, minOK=0, rawBad=0, minBad=0, dropped=0, minEqG=0)
for seed in range(1, seeds + 1):
    rng = random.Random(seed)
    for _ in range(per):
        n = rng.choice([5, 6, 7])
        Q, rv = G.randomInstance(rng, n)
        if not rv:
            continue
        for i in sorted(Q.nodes):
            if not ap.arrowsOutOf(Q, i):
                continue
            r = G.check(Q, rv, i)
            if not r["r7"] > r["g"]:
                continue
            t["excess"] += 1
            if not pr.isMutable(Q, rv, i):
                t["refused"] += 1
                continue
            nq, nr = pr.mutateAtVertex(Q, rv, i)
            if not nx.is_directed_acyclic_graph(nq):
                t["cyc"] += 1
                continue
            t["admitted"] += 1
            nr = [dict(x) for x in nr if x]
            mn = pr._minimalGenerators(nq, [dict(x) for x in nr])
            z = lambda rel: not any(any(row) for row in pr.cartanDiscrepancy(Q, rv, i, nq, rel))
            ok1, ok2 = z(nr), z(mn)
            t["rawOK"] += ok1; t["rawBad"] += (not ok1)
            t["minOK"] += ok2; t["minBad"] += (not ok2)
            t["dropped"] += sum(1 for x in nr if ap.source(x) == i) - sum(1 for x in mn if ap.source(x) == i)
            t["minEqG"] += sum(1 for x in mn if ap.source(x) == i) == r["g"]
print(t)
