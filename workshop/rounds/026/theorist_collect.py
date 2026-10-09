"""Round 026 (theorist): collect out-degree-2 mutable rows where both out-arrows carry a relation ending through v, with data for analysis.
Usage: theorist_collect.py n budget_sec class out.json   (walk as in rounds/025/theorist_control.py)"""
import sys, time, json
sys.path.insert(0, '.'); MYARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec'))
sys.argv = MYARGV
n, budget, cls, out = int(MYARGV[1]), float(MYARGV[2]), int(MYARGV[3]), MYARGV[4]
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; recs = []
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
nalg = 0
while frontier:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget: frontier.clear(); break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            outs = ap.arrowsOutOf(alg.quiver, v)
            if len(outs) == 2:
                rels = procedure.relationsFrom(alg)
                r = row(alg, v)
                recs.append(dict(v=v, arrows=sorted(alg.quiver.edges()), rels=[[list(map(list, q)) if False else [list(p) for p in rel]] for rel in alg.rels] if False else [[list(p) for p in rel] for rel in alg.rels],
                                 outs=[b[1] for b in outs], kerdim=r['kerdim'], mono=r['mono'], longsq=r['longsq']))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
json.dump(dict(n=n, cls=cls, algebras=nalg, recs=recs), open(out, 'w'))
print('algebras', nalg, 'rows', len(recs), '%.0fs' % (time.time() - t0))
