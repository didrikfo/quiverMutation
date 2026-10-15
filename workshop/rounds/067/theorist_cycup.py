"""Powered version of theorist_cycoff.py (round 067, theorist): take a round-065 acyclic instance with relations
(randomInstance, n = 5..7) and ATTACH an oriented 2-cycle c1 <-> c2 (labels n+1, n+2) that is not reachable from any
old vertex: extra arrows c -> u into old vertices (so the cycle is upstream of i) and/or extra relations starting at the
cycle and passing through old vertices (c -> u -> w) or through i.  Then for every old i with an out-arrow compare
K_opp = red(Z_1 e_i) with the repo's raw step-7 kernel and g with gd (same code as theorist_g2exact.py `check`).
Run: timeout 10m .venv/bin/python -W ignore workshop/rounds/067/theorist_cycup.py [seeds] [per_seed]
"""
import sys, random, importlib.util
sys.path.insert(0, ".")
from quivermutation import procedure as pr, arrowPaths as ap
spec = importlib.util.spec_from_file_location("g2", "workshop/rounds/065/theorist_g2exact.py")
g2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(g2)

if __name__ == "__main__":
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 250
    T = dict(pairs=0, bad=0, kpairs=0, kmulti=0, g=0, gd=0, gneqgd=0, exc=0, thruI=0, par=0, minEq=0, minNe=0)
    for seed in range(1, seeds + 1):
        rng = random.Random(seed)
        for _ in range(per):
            n = rng.choice([5, 6, 7])
            Q, rv = g2.randomInstance(rng, n)
            if not rv: continue
            c1, c2 = n + 1, n + 2
            Q.add_nodes_from([c1, c2]); Q.add_edge(c1, c2, key=0); Q.add_edge(c2, c1, key=0)
            for c in (c1, c2):
                for u in rng.sample(range(1, n + 1), rng.randint(0, 2)):
                    Q.add_edge(c, u, key=0)
            # relations starting on the cycle and running through old vertices
            for _ in range(rng.randint(0, 3)):
                ps = [p for p in ap.allPathsBetween(Q, rng.choice([c1, c2]), rng.randint(1, n)) if len(p) >= 2]
                if ps: rv.append(ap.combination([(rng.choice(ps), 1)]))
            par = ap.hasParallelArrows(Q)
            # relations that start on the new cycle: the set of vertices their paths visit after the start
            cycrels = [set(a[1] for t in rel for a in t) for rel in rv
                       if any(t and t[0][0] in (c1, c2) for t in rel)]
            for i in range(1, n + 1):
                if not ap.arrowsOutOf(Q, i): continue
                try:
                    r = g2.check(Q, rv, i)
                except Exception as e:
                    T["exc"] += 1; print("EXC", seed, i, repr(e)); continue
                T["pairs"] += 1; T["par"] += par
                T["thruI"] += any(i in vs for vs in cycrels)
                T["bad"] += (not r["ok"]); T["kpairs"] += r["kdim"] > 0; T["kmulti"] += r.get("kmulti", 0) > 0
                T["g"] += r["g"]; T["gd"] += r["gd"]; T["gneqgd"] += r["g"] != r["gd"]
                T["minEq"] += r["fromIstar"] == r["g"]; T["minNe"] += r["fromIstar"] != r["g"]
                if not r["ok"]: print("DIFFER", seed, i, sorted(Q.edges(keys=True)), rv)
    print(T)
    print("ALL EQUAL, g = gd" if T["bad"] == 0 and T["gneqgd"] == 0 and T["exc"] == 0 else "SOME DIFFER")
