"""Round 115 (theorist): HH^0,1,2 on ALL J != 0 steps recorded by the n = 7 class-1 key-guarded walk (checkpoint /tmp/skeptic_c1_13_ck1.pkl): kept steps
(child has the parent's key) and refused steps (child key differs: parent and child are then NOT derived equivalent, Coxeter polynomial / Cartan differs)."""
import sys, pickle, time
sys.path.insert(0, '.'); sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
exec(open('workshop/rounds/115/theorist_hh.py').read().split("if __name__")[0])
seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
recs = {}
for r in list(byCart.values()) + list(byKey.values()): recs[id(r)] = r
recs = list(recs.values()); print('J!=0 records', len(recs), 'nexp', nexp)
def H(alg):
    return hh(sorted(alg.quiver.nodes), sorted(alg.quiver.edges(keys=True)), procedure.relationsFrom(alg), top=2, maxdim=4000)
cnt = {}
for i, r in enumerate(recs):
    try:
        hp = H(r['alg']); hc = H(r['ch'])
    except Exception as ex:
        print(i, 'error', repr(ex)[:80]); continue
    if hp[1] is None or hc[1] is None: cnt['skipped(cap)'] = cnt.get('skipped(cap)', 0) + 1; continue
    k = ('kept' if r['kept'] else 'refused', 'HH equal' if hp[1] == hc[1] else 'HH differ')
    cnt[k] = cnt.get(k, 0) + 1
    if hp[1] != hc[1] or hp[1] != [1, 0, 0]: print(i, 'kept' if r['kept'] else 'refused', 'v', r['v'], 'parent', hp[1], 'child', hc[1], 'D0', r['D0'] if 'D0' in r else '?')
for k, v in sorted(cnt.items(), key=str): print(k, v)
