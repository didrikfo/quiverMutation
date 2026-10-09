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
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; hits = []; tab = Counter(); parentsKeys = set()
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
            if r['kerdim'] and not r['longsq']:
                hits.append((alg, v, r)); tab[(r['outdeg'],)] += 1
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('n', n, 'class', cls, 'algebras', len(seen), '%.0fs' % (time.time() - t0), 'rejections', len(hits), dict(tab))
cox = Counter(); dist = set()
for alg, v, r in hits:
    ck = search._coxeterKeyOrNone(alg); cox[ck == base] += 1
    dist.add(fingerprint.canonicalKey(alg))
print('parent Coxeter key == class base:', dict(cox), 'distinct parents (canonicalKey):', len(dist), 'distinct (parent,v) rows', len(hits))
shape = Counter()
for alg, v, r in hits:
    outs = ap.arrowsOutOf(alg.quiver, v)
    rels = alg.rels
    # per out arrow: does a relation of alg.rels end in that arrow (through v)?
    per = [sum(1 for rel in rels if any(len(q) >= 2 and q[-1] == b[1] and q[-2] == v for q in rel)) for b in outs]
    shape[(r['outdeg'], tuple(per), tuple(sorted(len(q) for rel in rels for q in rel)), r['mono'], len(alg.vertices()))] += 1
for k, x in shape.most_common(): print(k, x)
for alg, v, r in hits[:3]:
    print('SAMPLE v=', v, 'arrows', sorted(alg.quiver.edges()), 'rels', alg.rels, 'relsFrom', procedure.relationsFrom(alg), r)
# presentation test: is some relation of alg.rels redundant (in the ideal of the others)?  and does the rejection survive removing redundant ones?
def redundant(alg):
    out = []
    for j, rel in enumerate(alg.rels):
        others = [r for k, r in enumerate(alg.rels) if k != j]
        rr = procedure.relationsFrom(alg)
        oth = [x for k, x in enumerate(rr) if k != j]
        try:
            ok = ap.isInIdeal(alg.quiver, oth, rr[j])
        except Exception as e:
            ok = None
        if ok: out.append(j)
    return out
nred = Counter()
for alg, v, r in hits:
    red = redundant(alg)
    if red:
        keep = [x for k, x in enumerate(procedure.relationsFrom(alg)) if k not in red]
        # (removing several at once could change I; verify by kerdim with the pruned list only if I unchanged)
        same = all(ap.isInIdeal(alg.quiver, keep, x) for x in procedure.relationsFrom(alg))
        kd = kerdim(alg, v, keep)[0] if same else None
        nred[('redundant', len(red), 'pruned still generates I', same, 'kerdim after prune', kd)] += 1
    else: nred['minimal'] += 1
for k, x in nred.items(): print('PRESENTATION', k, x)
