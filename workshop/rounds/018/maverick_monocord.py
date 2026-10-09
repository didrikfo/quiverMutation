"""H-017 / T6 (round 018, maverick): can a MONOMIAL algebra with a cord (n arrows, connected, one undirected cycle,
no parallel arrows, acyclic quiver) share its Coxeter polynomial with ANY linear Nakayama algebra of n vertices?
Necessary condition for being in an LNA's mutation orbit. Enumerates every unicyclic acyclic quiver up to isomorphism
and every monomial ideal (antichain of paths of length >= 2).  usage: maverick_monocord.py N   (repository root)"""
import sys, os, itertools, collections
sys.path.insert(0, os.getcwd())
import networkx as nx
from quivermutation import coxeterTables as ct, invariants
n = int(sys.argv[1])
pairs = list(itertools.combinations(range(n), 2))
def quivers():
    buckets = collections.defaultdict(list)
    for es in itertools.combinations(pairs, n):
        g = nx.Graph(es)
        if g.number_of_nodes() != n or not nx.is_connected(g): continue
        for ori in itertools.product([0, 1], repeat=n):
            arr = tuple((a, b) if o == 0 else (b, a) for (a, b), o in zip(es, ori))
            d = nx.DiGraph(arr)
            if not nx.is_directed_acyclic_graph(d): continue
            h = nx.weisfeiler_lehman_graph_hash(d, iterations=3)
            if any(nx.is_isomorphic(d, e) for e, _ in buckets[h]): continue
            buckets[h].append((d, arr))
    return [arr for b in buckets.values() for _, arr in b]
def paths_of(arr):
    out = {v: [] for v in range(n)}
    for a, b in arr: out[a].append(b)
    res = []
    def go(p):
        if len(p) >= 3: res.append(tuple(p))
        for w in out[p[-1]]: go(p + [w])
    for a, b in arr: go([a, b])
    return res
def contains(p, q):  # q contiguous in p
    return any(p[i:i + len(q)] == q for i in range(len(p) - len(q) + 1))
lna = ct.lnaKeyIndex(n)
qs = quivers(); print("n", n, "unicyclic acyclic quivers", len(qs), "LNA polynomials", len(lna), flush=True)
tot = 0; hits = []; polys = set()
for arr in qs:
    P = paths_of(arr)  # vertex-sequence paths of length >= 2 arrows
    # antichain subsets under 'contains'
    def rec(i, chosen):
        if i == len(P): yield list(chosen); return
        yield from rec(i + 1, chosen)
        if all(not contains(P[i], c) and not contains(c, P[i]) for c in chosen):
            chosen.append(P[i]); yield from rec(i + 1, chosen); chosen.pop()
    allp = []
    out = {v: [] for v in range(n)}
    for a, b in arr: out[a].append(b)
    def gen(p):
        allp.append(tuple(p))
        for w in out[p[-1]]: gen(p + [w])
    for v in range(n): gen([v])
    for rel in rec(0, []):
        M = [[0] * n for _ in range(n)]
        for p in allp:
            if not any(contains(p, r) for r in rel): M[p[-1]][p[0]] += 1
        tot += 1
        key = invariants.coxeterCoefficients(M); polys.add(key)
        if key in lna: hits.append((arr, rel, lna[key]))
print("monomial cord algebras", tot, "distinct polynomials", len(polys), "sharing a polynomial with an LNA:", len(hits))
for arr, rel, l in hits[:12]: print(" arrows", arr, "rels", rel, "LNAs", ["".join(map(str, x)) for x in l][:4])
