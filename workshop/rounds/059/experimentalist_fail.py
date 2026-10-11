"""Independent n=6 check of the failure mode: for every LNA of length n and every right mutation (v in 1..n)
that is possible but NOT an LNA, classify the result quiver (nodes, edges, parallel arrows, cycles, relation shapes) by
the reason class of the rule (r_v>=3 / arrow v left-covered / r_v=0 and v<n-1 / v=n).  Usage: fail.py N"""
import sys, collections, networkx as nx
from quivermutation import lnaMoves as lm, mutation, nakayama
n = int(sys.argv[1])
rows = __import__('json').load(open(f"workshop/rounds/059/experimentalist_which_{n}.json"))
def cls(r, v):
    if v == n: return 'v=n'
    if any(r[j-1] and j + r[j-1] - 1 >= v for j in range(1, v) if j <= n-2): return 'left-covered'
    if v <= n-2 and r[v-1] >= 3: return 'r_v>=3'
    if v <= n-2 and r[v-1] == 0: return 'r_v=0'
    return 'other'
out = collections.Counter(); ok = 0
for r, res in rows:
    alg = nakayama.LinearNakayamaAlgebra(n, r)
    good = {v for v, _ in res.get('R', [])}
    for v in range(1, n + 1):
        if not lm._quiet(mutation.mutationIsPossibleAtVertex, alg, v): out[(cls(r, v), 'impossible')] += 1; continue
        nx_ = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [v])
        if v in good: ok += 1; continue
        q = nx_.quiver
        multi = q.number_of_edges() - len(set(q.edges()))  # parallel arrows
        cyc = bool(list(nx.simple_cycles(q))) if q.number_of_nodes() < 12 else None
        nrel = len(nx_.rels); sing = all(len(x) == 1 for x in nx_.rels)
        tag = f"nodes={len(q.nodes)} edges={q.number_of_edges()} parallel={multi} cycle={cyc} rels={nrel} allMonomial={sing}"
        out[(cls(r, v), tag)] += 1
print("n", n, "possible & LNA:", ok)
for k, c in sorted(out.items()): print(c, k)
