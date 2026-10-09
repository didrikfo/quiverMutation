"""T10(ii) round 049: read a checkpoint written by experimentalist_deepreplay.py (round 049), contract the seed LNAs
(+duals) to one node s, and ask of each failing key-keeping edge e=(u,w): (A) is e in the biconnected block of s
(then e lies on a simple path between two seeds), (B) walk length d_{G-e}(s,u)+1+d_{G-e}(s,w) against L*, the same
quantity minimised over the merge edges (child already seen). Explored graph only: unexpanded nodes are leaves.
usage: experimentalist_blocks.py CKPT"""
import sys, pickle
from collections import Counter, deque
import networkx as nx
sys.path.insert(0, '.')
S = pickle.load(open(sys.argv[1], 'rb')); E = S['edges']
firstChild, firstPar = {}, {}
for i, (p, c, f) in enumerate(E):
    firstPar.setdefault(p, i); firstChild.setdefault(c, i)
seeds = {p for p in firstPar if p not in firstChild or firstPar[p] < firstChild[p]}
print('edges', len(E), 'nodes', len(set(x for e in E for x in e[:2])), 'seeds', len(seeds), 'expanded', S['expanded'], 'depth', S['d'], 'pos', S['pos'])
def node(x): return 's' if x in seeds else x
G = nx.Graph(); mult = Counter(); fails = []
for p, c, f in E:
    a, b = node(p), node(c)
    if a == b: continue
    k = frozenset((a, b)); mult[k] += 1; G.add_edge(a, b)
    if f: fails.append((a, b))
print('contracted graph', G.number_of_nodes(), G.number_of_edges(), 'failing edges', len(fails))
blocks = [set(b) for b in nx.biconnected_components(G) if 's' in b]
def inblock(a, b): return any(a in B and b in B for B in blocks)
def bfs(skip):
    d = {'s': 0}; q = deque(['s'])
    while q:
        u = q.popleft()
        for w in G[u]:
            if frozenset((u, w)) == skip or w in d: continue
            d[w] = d[u] + 1; q.append(w)
    return d
d0 = bfs(None)
# merge edges: edges between nodes where a (non-tree) cycle exists; L* over edges whose removal keeps both ends connected
Lmin = None; Lcnt = Counter()
for a, b in G.edges():
    if a == 's' or b == 's':
        if mult[frozenset((a, b))] == 1 and (a == 's' and b == 's'): continue
    if abs(d0[a] - d0[b]) == 0 or True:
        pass
# merge edge = non-tree BFS edge: |d(a)-d(b)|<=0 (same level) or an extra parent
par = Counter()
for a, b in G.edges():
    if d0[a] + 1 == d0[b]: par[b] += 1
    elif d0[b] + 1 == d0[a]: par[a] += 1
merge = [(a, b) for a, b in G.edges() if d0[a] == d0[b] or (d0[a] + 1 == d0[b] and par[b] > 1) or (d0[b] + 1 == d0[a] and par[a] > 1)]
Ls = [d0[a] + 1 + d0[b] for a, b in merge]
print('merge edges (same level or second parent):', len(merge), 'L* (shortest merge walk length) =', min(Ls), 'histogram', sorted(Counter(Ls).items())[:8])
print('failing edges:')
tot = Counter()
for a, b in fails:
    k = frozenset((a, b)); d = bfs(k)
    L = d[a] + 1 + d[b] if a in d and b in d else None
    ib = inblock(a, b); m = mult[k]
    tot[(ib, m > 1, L is not None)] += 1
    print('  depth(u,w) %s,%s  in-degree(w)=%d deg(w)=%d  inBlockOfSeeds=%s  mult=%d  walkLen=%s' % (d0[a], d0[b], par[b], G.degree(b), ib, m, L))
print('summary (inBlock, multiple, connectedWithout):', dict(tot))
