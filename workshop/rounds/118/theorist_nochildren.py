"""Round 118 (theorist): E-305 open point -- are the 2 NO children (child HH^1 = 1) of the refused J != 0 steps tree quivers?
Needs /tmp/skeptic_c1_13_ck1.pkl (rebuild: `timeout 10m .venv/bin/python -u workshop/rounds/069/skeptic_c1_13.py 1 540 23000`, twice).
Prints, for every refused step: undirected cycle rank of parent and child quiver (|E| - |V| + #components); for those whose child HH differs from the parent,
the child's edges and relations.  Tree quiver <=> cycle rank 0 (and then HH^m = 0, m >= 1, is classical -- which is the check on the HH code)."""
import sys, pickle
sys.path.insert(0, '.'); sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
exec(open('workshop/rounds/115/theorist_hh.py').read().split("if __name__")[0])
seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
recs = {}
for r in list(byCart.values()) + list(byKey.values()): recs[id(r)] = r
recs = list(recs.values()); print('J!=0 records', len(recs), 'nexp', nexp)
def H(alg): return hh(sorted(alg.quiver.nodes), sorted(alg.quiver.edges(keys=True)), procedure.relationsFrom(alg), top=2, maxdim=4000)
def rank(alg):
    V = sorted(alg.quiver.nodes); E = sorted(alg.quiver.edges(keys=True)); par = {v: v for v in V}
    def f(x):
        while par[x] != x: x = par[x]
        return x
    for s, t, k in E: par[f(s)] = f(t)
    return len(E) - len(V) + len({f(v) for v in V})
from collections import Counter
cnt = Counter()
for i, r in enumerate(recs):
    if r['kept']: continue
    hp, hc = H(r['alg']), H(r['ch'])
    rp, rc = rank(r['alg']), rank(r['ch'])
    cnt[('child HH', tuple(hc[1]) if hc[1] else None, 'parent cycle rank', rp, 'child cycle rank', rc)] += 1
    if hp[1] != hc[1]:
        print('NO child', i, 'v', r['v'], 'parent HH', hp[1], 'child HH', hc[1], 'parent cycle rank', rp, 'child cycle rank', rc)
        print('  parent edges', sorted(r['alg'].quiver.edges(keys=True)), '\n  parent rels', procedure.relationsFrom(r['alg']))
        print('  child  edges', sorted(r['ch'].quiver.edges(keys=True)), '\n  child  rels', procedure.relationsFrom(r['ch']))
for k, v in sorted(cnt.items(), key=str): print(v, k)
