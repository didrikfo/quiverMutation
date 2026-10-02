"""Round 025 (skeptic): collect the n=8 class-0 gate-admitted rejections with out-degree 2 and no long square (round 023 scholar_longsquare walk),
print Coxeter-class membership of the parent, structure (shared relations, minimality, out-arrow relations) and a sample.
Usage: skeptic_n8.py n budget_sec [class]"""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); MYARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec'))
sys.argv = MYARGV
n, budget = int(MYARGV[1]), float(MYARGV[2]); cls = int(MYARGV[3]) if len(MYARGV) > 3 else 0
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; hits = []; ctl=Counter(); tab = Counter(); parentsKeys = set()
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
while frontier:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget: frontier.clear(); break
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            r = row(alg, v)
            if r['outdeg']==2:
                outs=ap.arrowsOutOf(alg.quiver, v)
                per=[sum(1 for rel in alg.rels if any(len(q) >= 2 and q[-1] == b[1] and q[-2] == v for q in rel)) for b in outs]
                ctl[(bool(r['kerdim']), all(per))]+=1
            if r['kerdim'] and not r['longsq']:
                hits.append((alg, v, r)); tab[(r['outdeg'],)] += 1
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('CONTROL (kerdim>0, both arrows carry a relation) among gate-or-not outdeg2 rows:',dict(ctl))
