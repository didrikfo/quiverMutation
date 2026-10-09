"""Round 038 (toolsmith): is the key-preserving n = 6 class infinite? Per BFS level, max arrow count and max parallel-arrow multiplicity.
Usage: toolsmith_n6growth.py CLASSIDX LEVELS   (class 0 = the 2-LNA key class holding the E-134 key coincidences)"""
import sys, collections
a = sys.argv[1:]; sys.argv = ['x']
src = open('workshop/rounds/038/toolsmith_n6close.py').read().split("hits = json.load")[0]
sys.argv = ['x', '--plan', '1']; exec(compile(src, 'c', 'exec'))
idx, L = int(a[0]), int(a[1]); key = order[idx]
seen = set(); fr = []
def add(x):
    ck = fingerprint.canonicalKey(x)
    if ck in seen: return
    seen.add(ck); fr.append(x)
for x in classes[key]: add(x)
for lvl in range(1, L + 1):
    cur, fr[:] = fr[:], []; mx = collections.Counter(); par = 0
    for alg in cur:
        E = list(alg.quiver.edges()) if not hasattr(alg.quiver, 'edges') or True else []
        c = collections.Counter(E); mx['arrows'] = max(mx['arrows'], len(E)); par = max(par, max(c.values()) if c else 0)
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == key: add(ch)
    print('level', lvl, 'expanded', len(cur), 'max arrows', mx['arrows'], 'max parallel mult', par, 'seen', len(seen), flush=True)
