"""Round 091 response: per-depth flagged counts of the 2400230 ball scan, and key(2400230) vs base key.
usage (repo root): experimentalist_jscan2.py STATE.pkl [MAXDEPTH]"""
import sys, pickle, collections
ST = sys.argv[1]; MAXD = int(sys.argv[2]) if len(sys.argv) > 2 else 8
sys.argv = ['x', 'ctl', 'x']
src = open('workshop/rounds/087/toolsmith_d6scan.py').read().split("if MODE == 'build'")[0]
sys.path.insert(0, '.'); exec(compile(src, 's', 'exec'))
for c in ('2400230', '3060000', '3304000'):
    print('key', c, search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName(c)) == base)
S = pickle.load(open(ST, 'rb')); seen = S['seen']
nodes = collections.Counter(); flagged = collections.Counter(); events = collections.Counter(); kept = collections.Counter()
for k, (x, pk, mv, d) in seen.items():
    if d > MAXD: continue
    nodes[d] += 1; c = jcands(x)
    if c: flagged[d] += 1
    for kind, v, st in c:
        events[d] += 1
        if st == 'kept': kept[d] += 1
for d in range(MAXD + 1): print('depth', d, 'nodes', nodes[d], 'parents_flagged', flagged[d], 'steps', events[d], 'kept', kept[d])
print('total flagged parents', sum(flagged.values()), 'steps', sum(events.values()))
