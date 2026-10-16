"""Round 069 (skeptic), T10/T5: rebuild the 13 class-1 E-147 steps (n = 7) independently of skeptic_c2.py.
Own walk (key-guarded BFS from the LNAs + duals of key class CLS, canonicalKey dedup, as E-151/E-147), own det-polynomial
D(x) = det(xC_B + C_B^T) - det(xC_A + C_A^T) (sympy, exact), own tilting test (rounds/050/skeptic_tilt.py: Hom_K(T,T[m]), m=-1,0,1,
not tiltingPlus/perI/gate) on every J != 0 step.  Records: all gate-admitted legal vertex-preserving J != 0 steps, distinct by
(Cartan A, Cartan B, v) [E-147 convention] and by (canonicalKey(parent), v) [E-151 convention]; D == 0; key kept; tiltingPlus; parallel
arrows in the parent; membership in E-151's set (key kept and tiltingPlus false).
Usage: skeptic_c1_13.py CLS DEADLINE_SEC MAXEXP       (run from repo root; every claim is under the J = 0 premise)"""
import sys, time, pickle
from collections import Counter
import numpy as np, sympy
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
CLS, DEAD, MAXEXP = int(_a[1]), float(_a[2]), int(_a[3]); n = 7
from fractions import Fraction
def idet(M):
    """exact integer determinant (Fraction elimination)."""
    M = [[Fraction(int(x)) for x in r] for r in M]; m = len(M); d = Fraction(1)
    for c in range(m):
        p = next((r for r in range(c, m) if M[r][c] != 0), None)
        if p is None: return 0
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, m):
            f = M[r][c] / M[c][c]
            if f: M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return int(d)
_dp = {}
def detpoly(C):
    """det(xC + C^T) is a polynomial of degree <= 7 in x: its values at 9 integers determine it exactly."""
    k = C.tobytes()
    if k not in _dp: _dp[k] = tuple(idet(x * C + C.T) for x in range(-4, 5))
    return _dp[k]
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[CLS]
print('class sizes', [len(classes[k]) for k in order], 'using', CLS, 'seeds', len(classes[base]), flush=True)
from collections import deque
CK = '/tmp/skeptic_c1_13_ck%d.pkl' % CLS
import os
if os.path.exists(CK):
    seen, queue, byCart, byKey, nexp, snaps = pickle.load(open(CK, 'rb')); print('resumed at nexp', nexp, flush=True)
else:
    seen = set(); queue = deque(); byCart = {}; byKey = {}; nexp = 0; snaps = {}
    for s in classes[base]:
        k = fingerprint.canonicalKey(s)
        if k is None or k not in seen:
            if k is not None: seen.add(k)
            queue.append((s, 0))
t0 = time.time(); stop = False
while queue:
    if nexp >= MAXEXP or time.time() - t0 > DEAD: stop = True; break
    alg, lv = queue.popleft(); nexp += 1
    if nexp in (2500, 5000, 10000, 15000, 20000): snaps[nexp] = (len(byCart), len(byKey))
    if list(nx.simple_cycles(alg.quiver)): continue
    V = sorted(alg.quiver.nodes); rels = procedure.relationsFrom(alg)
    for v in V:
        if not mutation.mutationIsPossibleAtVertex(alg, v): continue
        raw = mutation.quiverMutationAtVertex(alg, v)
        if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
        ch = reduction.reducePathAlgebra(raw)
        if sorted(ch.vertices()) != V: continue
        ck = search._coxeterKeyOrNone(ch)
        J = perI(alg, v)
        if J:
            CA = np.array(cartan(alg), dtype=int); CB = np.array(cartan(ch), dtype=int); vi = V.index(v)
            ckey = fingerprint.canonicalKey(alg); k1 = (CA.tobytes(), CB.tobytes(), vi); k2 = (ckey, v)
            if k1 not in byCart or k2 not in byKey:
                tp = bool(tiltingPlus(alg.quiver, rels, v))
                rec = dict(alg=alg, ch=ch, v=v, V=V, J=J, depth=lv, kept=(ck == base), tp=tp, nexp=nexp, CA=CA, CB=CB,
                           par=sorted({(a, b) for a, b, kk in alg.quiver.edges(keys=True) if alg.quiver.number_of_edges(a, b) > 1}))
                if k1 not in byCart: byCart[k1] = rec
                if k2 not in byKey: byKey[k2] = rec
        if ck == base:
            k = fingerprint.canonicalKey(ch)
            if k is None or k not in seen:
                if k is not None: seen.add(k)
                queue.append((ch, lv + 1))
pickle.dump((seen, queue, byCart, byKey, nexp, snaps), open(CK, 'wb'))
print('segment end: nexp', nexp, 'queue', len(queue), 'distinct Cartan', len(byCart), 'distinct key', len(byKey), '%.0fs' % (time.time() - t0), flush=True)
if queue and nexp < MAXEXP: sys.exit(0)
stop = bool(queue)
print('walk: expansions', nexp, 'stopped by cap', stop, 'snapshots (distinct Cartan, distinct key, sec)', snaps, flush=True)
for r in list(byCart.values()) + list(byKey.values()):
    if 'D0' not in r:
        dA, dB = detpoly(r['CA']), detpoly(r['CB']); r['D0'] = (dA == dB)
# E-147 convention
S147 = list(byCart.values()); Z = [r for r in S147 if r['D0']]
print('E-147 convention: distinct J!=0 steps', len(S147), 'D==0:', len(Z), 'of which key kept', sum(r['kept'] for r in Z), 'tiltingPlus false', sum(not r['tp'] for r in Z), 'parent has parallel arrows', sum(bool(r['par']) for r in Z))
S151 = [r for r in byKey.values() if r['kept'] and not r['tp']]
print('E-151 convention: distinct (parent,v) key kept & tiltingPlus false:', len(S151), ' D==0 among them:', sum(r['D0'] for r in S151), ' parallel parents:', sum(bool(r['par']) for r in S151))
print('all distinct (parent,v) J!=0:', len(byKey), ' key kept:', sum(r['kept'] for r in byKey.values()), ' D==0 & kept:', sum(r['D0'] and r['kept'] for r in byKey.values()))
k151 = {id(r) for r in S151}
# independent tilting test on E-147's 13 and on E-151's 16
def indep(r):
    V, res = stepTest(r['alg'], r['v'], -1)
    m1 = int(np.array(res[-1]).sum()); p1 = int(np.array(res[1]).sum()); return m1, p1
out = []
for tag, S in (('E147-D0', Z), ('E151-set', S151)):
    for i, r in enumerate(S):
        m1, p1 = indep(r)
        inE151 = any(r['alg'] is q['alg'] and r['v'] == q['v'] or (fingerprint.canonicalKey(r['alg']) == fingerprint.canonicalKey(q['alg']) and r['v'] == q['v']) for q in S151)
        row = (tag, i, 'depth', r['depth'], 'v', r['v'], 'J', dict(r['J']), 'tp', r['tp'], 'D0', r['D0'], 'kept', r['kept'], 'Hom(T,T[-1]) dim', m1, 'Hom(T,T[1]) dim', p1, 'parallel', r['par'], 'in E151', inE151, 'dim', int(r['CA'].sum()), int(r['CB'].sum()))
        print(*row, flush=True); out.append(row)
pickle.dump([dict(tag=t, v=r['v'], V=r['V'], pe=sorted(r['alg'].quiver.edges(keys=True)), pr=repr(procedure.relationsFrom(r['alg'])), ce=sorted(r['ch'].quiver.edges(keys=True)), cr=repr(procedure.relationsFrom(r['ch'])), D0=r['D0'], tp=r['tp'], kept=r['kept'], par=r['par']) for t, S in (('E147', Z), ('E151', S151)) for r in S], open('/tmp/skeptic_c1_13_cls%d.pkl' % CLS, 'wb'))
