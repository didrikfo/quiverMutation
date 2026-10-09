"""Round 026 (skeptic): for the n=8 class-0 walk (round-023 walk, same wall-clock cap), take every out-degree-2 gate-admitted row and
 (a) intrinsic control: dim J_beta (c e_i A e_v with c*beta=0) for each out-arrow and dim J = common kernel, per source vertex i;
 (b) extract a kernel element x (normal form in e_i A e_v) for each rejecting row (J != 0, outdeg 2);
 (c) shorter presentation: replace a relation by itself minus a term lying in the ideal of the others (451=0 turns 4513=4573 into 4573=0),
     recompute kerdim and the 'each out-arrow carries a relation' shape on the shortened presentation.
Usage: skeptic_x.py n budget_sec [class]"""
import sys, time
from collections import Counter
from fractions import Fraction
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
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; rows_ = []
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
            if r['outdeg'] == 2: rows_.append((alg, v, r))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('walk: algebras', len(seen), 'outdeg-2 rows', len(rows_), '%.0fs' % (time.time() - t0), flush=True)

def maprows(quiver, rels, i, v, P, outs):
    out = []
    for q in P:
        rw = {}
        for b in outs:
            for kk, x in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(quiver, rels, i, b[1])).items(): rw[(b, kk)] = x
        out.append(rw)
    return out
def jdim(quiver, rels, i, v, outs):
    P = ap.allPathsBetween(quiver, i, v)
    if not P: return 0, P
    dimV = len(P) - len(ap.idealBasis(quiver, rels, i, v))
    return dimV - rank(maprows(quiver, rels, i, v, P, outs)), P
def nullvec_residues(quiver, rels, i, v, outs):
    """normal forms (mod I) of null vectors of p -> (p beta)_beta on paths i~>v"""
    P = ap.allPathsBetween(quiver, i, v); rws = maprows(quiver, rels, i, v, P, outs); piv = {}; res = []
    for idx, r in enumerate(rws):
        r = {k: Fraction(x) for k, x in r.items()}; c = {idx: Fraction(1)}
        while r:
            h = min(r)
            if h in piv:
                f = r[h]; pr, pc = piv[h]
                r = {k: r.get(k, 0) - f * pr.get(k, 0) for k in set(r) | set(pr)}; r = {k: x for k, x in r.items() if x != 0}
                c = {k: c.get(k, 0) - f * pc.get(k, 0) for k in set(c) | set(pc)}; c = {k: x for k, x in c.items() if x != 0}
            else:
                f = r[h]; piv[h] = ({k: x / f for k, x in r.items()}, {k: x / f for k, x in c.items()}); break
        else:
            x = ap.reduceAgainstPivots(ap.combination({P[k]: x for k, x in c.items()}), ap.idealBasis(quiver, rels, i, v))
            if x: res.append(x)
    return res
def reducePresentation(quiver, rels):
    rels = [dict(r) for r in rels]; changed = True; nred = 0
    while changed:
        changed = False
        for idx, r in enumerate(rels):
            others = rels[:idx] + rels[idx + 1:]
            for pth in list(r):
                if len(r) == 1 and not others: continue
                if ap.isInIdeal(quiver, others, ap.combination([pth])):
                    if len(r) == 1: rels = others
                    else: rels[idx] = {k: x for k, x in r.items() if k != pth}
                    changed = True; nred += 1; break
            if changed: break
    return rels, nred
def shape(rels, outs):  # every out-arrow is the last arrow of a path in some relation
    return all(any(any(len(p) >= 1 and p[-1] == b for p in r) for r in rels) for b in outs)
def vl(p): return ''.join(str(a[0]) for a in p) + (str(p[-1][1]) if p else '')

ctl = Counter(); ctl2 = Counter(); hits = []; shp = Counter(); pers = Counter(); xs = []
for alg, v, r in rows_:
    q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(q, v)
    rej = bool(r['kerdim']) and not r['longsq']
    jb = [0, 0]; samei = False; jall = 0
    for i in q.nodes:
        if i == v: continue
        d = [jdim(q, rels, i, v, [b])[0] for b in outs]
        if d[0] > 0: jb[0] += 1
        if d[1] > 0: jb[1] += 1
        if d[0] > 0 and d[1] > 0: samei = True
    both_any = jb[0] > 0 and jb[1] > 0
    ctl[(bool(r['kerdim']), both_any, samei)] += 1
    rr, nred = reducePresentation(q, rels)
    sh0, sh1 = shape(rels, outs), shape(rr, outs)
    shp[(bool(r['kerdim']), sh0, sh1)] += 1
    if r['kerdim'] or sh0 != sh1:
        pass
    if rej:
        kd2 = sum(kerdim(alg, v, rr)[0:1])
        multi_before = sum(1 for x in rels if len(x) > 1); multi_after = sum(1 for x in rr if len(x) > 1)
        hits.append((alg, v, nred, len(rels), len(rr), kd2, multi_before, multi_after, sh0, sh1))
        for i in q.nodes:
            if i == v: continue
            for x in nullvec_residues(q, rels, i, v, outs):
                ok = all(not ap.reduceAgainstPivots(ap.combination({p + (b,): c for p, c in x.items()}), ap.idealBasis(q, rels, i, b[1])) for b in outs)
                xs.append((i, v, x, ok, (q, rels, outs))); break
print('CONTROL intrinsic: (J!=0, J_beta1!=0 and J_beta2!=0 (any i), same i):', dict(ctl))
print('SHAPE (J!=0, shape orig presentation, shape shortened):', dict(shp))
print('rejecting rows', len(hits), 'with a reducible relation', sum(1 for h in hits if h[2] > 0), 'kerdim>0 persists on shortened', sum(1 for h in hits if h[5] > 0),
      'shape persists', sum(1 for h in hits if h[9]), 'monomial-ised relations (multi->fewer terms)', sum(1 for h in hits if h[7] < h[6]))
seenx = {}
cnt = Counter()
for (i, v, x, ok, _) in xs:
    lens = sorted({len(p) for p in x}); coefs = sorted(set(x.values()))
    cnt[(len(x), tuple(lens), tuple(str(c) for c in coefs), ok)] += 1
print('KERNEL ELEMENT x (one per rejecting row): (terms, path lengths, coefficients, x*beta in I for all beta):')
for k, c in sorted(cnt.items(), key=str): print(' ', k, c)
kill = Counter()
for (i, v, x, ok, (q, rels, outs)) in xs:
    per = []
    for b in outs:
        per.append(sum(1 for pth in x if not ap.reduceAgainstPivots(ap.combination([pth + (b,)]), ap.idealBasis(q, rels, i, b[1]))))
    kill[tuple(sorted(per))] += 1
print('TERMS of x killed singly by each out-arrow (sorted pair; 0 = only the sum is killed):', dict(kill))
print('x found for', len(xs), 'of', len(hits), 'rows')
for (i, v, x, ok, _) in xs[:6]:
    print('SAMPLE i=%d v=%d x =' % (i, v), ' + '.join('%s*%s' % (c, vl(p)) for p, c in x.items()))
for h in hits:
    if h[2] > 0:
        print('REDUCIBLE example: rels', h[3], '->', h[4], '; reducible steps', h[2]); break
