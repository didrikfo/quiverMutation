"""Round 026 (theorist): test the witness condition W on collected out-degree-2 rows.
W(alg,v): some relation R of alg (procedure.relationsFrom) all of whose terms end with the arrow b1 out of v and pass through v (q[-2]==v),
 x = R / b1 (prefixes with R's coefficients, a combination of paths i -> v) is not zero in A and x*b2 is zero in A for the other out-arrow b2.
Then x*b1 = 0 = x*b2, so x in J (kerdim > 0).  Compare W with kerdim>0 over all rows; split by the shape of x.
Usage: theorist_witness.py file.json"""
import sys, json
from collections import Counter
sys.path.insert(0, '.'); _a = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec')); sys.argv = _a
D = json.load(open(sys.argv[1])); tab = Counter(); mism = []
for r in D['recs']:
    nodes = sorted({x for e in r['arrows'] for x in e})
    A = build([tuple(e) for e in r['arrows']], r['rels'], nodes); v = r['v']; Q = A.quiver
    try: rels = procedure.relationsFrom(A)
    except ValueError: tab['parallel-skipped', bool(r['kerdim'])] += 1; continue
    outs = ap.arrowsOutOf(Q, v)
    kinds = []; 
    for R in rels:
        if len(R) < 2: continue
        for b1 in outs:
            if not all(len(q) >= 2 and q[-1] == b1 and q[-2][1] == v for q in R): continue
            b2 = [b for b in outs if b != b1][0]
            x = ap.combination({q[:-1]: c for q, c in R.items()})
            if ap.isInIdeal(Q, rels, x): continue
            xb2 = ap.combination({q[:-1] + (b2,): c for q, c in R.items()})
            if ap.isInIdeal(Q, rels, xb2):
                zero = all(ap.isInIdeal(Q, rels, ap.combination([q[:-1] + (b2,)])) for q in R)
                kinds.append((len(R), min(len(q) for q in R), 'termwise-zero' if zero else 'cancels'))
    W = bool(kinds); K = bool(r['kerdim']); tab[(W, K)] += 1
    if W != K: mism.append((r['v'], r['rels'], r['kerdim']))
    if K: tab[('rejshape',) + tuple(sorted(set(kinds)))] += 1
print(dict(tab)); print('mismatch', len(mism))
for m in mism[:6]: print(m)
