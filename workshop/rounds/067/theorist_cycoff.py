"""g = gd and Lemma L with an oriented cycle OFF the support of i (round 067, theorist).
Run from the repository root:  timeout 10m .venv/bin/python workshop/rounds/067/theorist_cycoff.py [seeds] [per_seed]

Instances: random multigraph on 1..n (forward edges, parallel arrows allowed) plus one or two BACK edges
b->a (a<b) that close oriented cycles; relations are random monomial/binomial combinations of
vertex-simple paths of length >= 2.  For each vertex i with an out-arrow we classify:
  'off'  : no vertex reachable from i (including i) lies on an oriented cycle  -> the descendant set D(i) is acyclic
  'on'   : some cycle meets D(i)  (cycle downstream of, or through, i)
and, for class 'off', compare (as in round 065) K_opp = red(Z_1 e_i) with the repo's raw step-7 kernel, and
g (=dim Z_1/(B_1 + Z_1 rad + rad Z_1)) with gd (=dim K/(sum K(w)a + red B_1)).  Class 'on' only counts g, gd with the
vertex-simple Opp (meaningless there: reported for contrast) and exceptions raised by the repo side.
"""
import sys, random, importlib.util
import networkx as nx
sys.path.insert(0, ".")
from quivermutation import procedure as pr, arrowPaths as ap
spec = importlib.util.spec_from_file_location("g2", "workshop/rounds/065/theorist_g2exact.py")
g2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(g2)


def inst(rng, n):
    Q = nx.MultiDiGraph(); Q.add_nodes_from(range(1, n + 1))
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            if rng.random() < 0.45:
                for k in range(2 if rng.random() < 0.2 else 1):
                    Q.add_edge(a, b, key=k)
    for _ in range(rng.choice([1, 1, 2])):
        a = rng.randint(1, n - 1); b = rng.randint(a + 1, n)
        Q.add_edge(b, a, key=0)
    rv = []
    for _ in range(rng.randint(1, 6)):
        a = rng.randint(1, n); b = rng.randint(1, n)
        ps = [p for p in ap.allPathsBetween(Q, a, b) if len(p) >= 2]   # vertex-simple, a != b or simple cycle
        if not ps: continue
        if len(ps) >= 2 and rng.random() < 0.4:
            p1, p2 = rng.sample(ps, 2)
            if ap.source(ap.combination([(p1, 1)])) != ap.source(ap.combination([(p2, 1)])): continue
            rv.append(ap.combination([(p1, 1), (p2, rng.choice([-1, -2, 2]))]))
        else:
            rv.append(ap.combination([(rng.choice(ps), 1)]))
    return Q, rv


def cyclicVertices(Q):
    S = set()
    for comp in nx.strongly_connected_components(nx.DiGraph(Q)):
        if len(comp) > 1 or any(Q.has_edge(v, v) for v in comp): S |= comp
    return S


if __name__ == "__main__":
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    T = dict(off=0, offbad=0, offg=0, offgd=0, offggtgd=0, offk=0, on=0, onexc=0, offexc=0, hascyc=0, upstream=0)
    for seed in range(1, seeds + 1):
        rng = random.Random(seed)
        for _ in range(per):
            n = rng.choice([5, 6, 7])
            Q, rv = inst(rng, n)
            if not rv: continue
            cyc = cyclicVertices(Q)
            if not cyc: continue
            T["hascyc"] += 1
            opp = None
            for i in sorted(Q.nodes):
                if not ap.arrowsOutOf(Q, i): continue
                D = nx.descendants(Q, i) | {i}
                if D & cyc:
                    T["on"] += 1
                    continue
                T["off"] += 1
                T["upstream"] += any(nx.has_path(Q, c, i) for c in cyc)
                try:
                    r = g2.check(Q, rv, i)
                except Exception as e:
                    T["offexc"] += 1; print("EXC off", seed, i, repr(e), sorted(Q.edges(keys=True)), rv); continue
                T["offbad"] += (not r["ok"]); T["offg"] += r["g"]; T["offgd"] += r["gd"]
                T["offggtgd"] += r["g"] != r["gd"]; T["offk"] += r["kdim"] > 0
                if not r["ok"]: print("DIFFER", seed, i, sorted(Q.edges(keys=True)), rv)
    print(T)
    print("ALL EQUAL, g = gd" if T["offbad"] == 0 and T["offggtgd"] == 0 and T["offexc"] == 0 else "SOME DIFFER")
