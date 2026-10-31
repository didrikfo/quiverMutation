"""Round 103 (theorist), T5: among small finite-dimensional monomial algebras with an oriented cycle, how often does the library's gate
(procedure.isMutable) admit a vertex i that lies on an oriented cycle (the only way a loop of degree >= 1 at i can occur in
Oppermann's model: a loop at i in the model forces a closed path through i in Q^(0))?  What do mutateAtVertex(..., allowCyclic=True) and the default do?
usage (repository root):  timeout 10m .venv/bin/python -W ignore workshop/rounds/103/theorist_cycsweep.py [seed] [samples]"""
import sys, random, collections, itertools
sys.path.insert(0, '.'); sys.path.insert(0, 'tests')
import networkx as nx
from quivermutation import procedure as pr, arrowPaths as ap
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1; N = int(sys.argv[2]) if len(sys.argv) > 2 else 400
rng = random.Random(seed)

def fd(arrows, rels, cap=14):
    """finite-dimensional?  (no nonzero path of length > cap): monomial, vertex-path relations given as arrow-index tuples."""
    zero = lambda p: any(any(p[k:k+len(r)] == r for k in range(len(p)-len(r)+1)) for r in rels)
    front = [(a,) for a in range(len(arrows)) if not zero((a,))]
    for L in range(2, cap + 1):
        front = [p + (a,) for p in front for a in range(len(arrows)) if arrows[p[-1]][1] == arrows[a][0] and not zero(p + (a,))]
        if not front: return True
    return False

stats = collections.Counter(); samples = []
tries = 0
while stats['fd_cyclic'] < N and tries < 200000:
    tries += 1
    n = rng.choice([3, 4])
    pairs = [(u, v) for u in range(1, n + 1) for v in range(1, n + 1) if u != v]
    arrows = rng.sample(pairs, rng.randint(n, min(len(pairs), n + 2)))
    G = nx.DiGraph(arrows)
    if nx.is_directed_acyclic_graph(G) or len(G.nodes) < n: continue
    paths = []   # arrow-index paths of length 2..4
    for L in (2, 3, 4):
        for p in itertools.product(range(len(arrows)), repeat=L):
            if all(arrows[p[k]][1] == arrows[p[k+1]][0] for k in range(L - 1)): paths.append(p)
    if not paths: continue
    rels = [p for p in rng.sample(paths, min(len(paths), rng.randint(1, 4)))]
    # minimal: drop relations containing another as a subword
    rels = [r for r in rels if not any(s != r and any(r[k:k+len(s)] == s for k in range(len(r)-len(s)+1)) for s in rels)]
    if not fd(arrows, rels): continue
    stats['fd_cyclic'] += 1
    Q = nx.MultiDiGraph(); Q.add_nodes_from(range(1, n + 1))
    for k, (t, h) in enumerate(arrows): Q.add_edge(t, h, key=0 if arrows.count((t, h)) == 1 else k)
    key = {k: (t, h, Q_k) for k, ((t, h), Q_k) in enumerate(zip(arrows, [0] * len(arrows)))}
    R = [ap.combination([(tuple(key[a] for a in r), 1)]) for r in rels]
    cyc = {v for c in nx.simple_cycles(nx.DiGraph(arrows)) for v in c}
    for i in range(1, n + 1):
        try: gate = pr.isMutable(Q, R, i)
        except Exception as e: stats['gate_raised'] += 1; continue
        if not gate: continue
        onc = i in cyc
        stats['gate_true_i_on_cycle' if onc else 'gate_true_i_off_cycle'] += 1
        for allow in (False, True):
            try: pr.mutateAtVertex(Q, R, i, allowCyclic=allow); out = 'returns'
            except Exception as e: out = type(e).__name__
            stats[('on_cycle' if onc else 'off_cycle', 'allowCyclic=%s' % allow, out)] += 1
        if onc and len(samples) < 3: samples.append((arrows, rels, i))
print('tries', tries); [print(' ', k, v) for k, v in sorted(stats.items(), key=repr)]
for s in samples: print('sample gate-admitted, i on a cycle:', s)
