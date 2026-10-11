"""Round 058 (skeptic), T10: equal-dims wrong-algebra control for the E-174 quiver-level End(T) test (symcheck2).
Needs /tmp/tsm/c1.pkl (rebuild: DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100, run twice).
Usage (repo root):  CLS=1 timeout 10m .venv/bin/python workshop/rounds/058/skeptic_wrongalg.py
For each failing c1 step whose child has one parallel pair P = (b0,b1) and 'line-type' relations (x then P, kill one linear combination of
x b0, x b1), re-choose the killed line of every such relation among 4 lines [1:0],[0:1],[1:1],[1:-1] (4^k variants).  For each variant:
 (i) dim of K Q/I' computed here by linear algebra (independent of the test) vs dim End(T) (= dims of the true child);
 (ii) ground-truth label-preserving iso with the true child: iff a Moebius map of P^1 carries the k killed points of c to those of c'
      (k <= 3 points: iff same coincidence pattern; this is an invariant argument, not a computation);
 (iii) symcheck2 verdict.  Reports the 2x2 table verdict x truth among EQUAL-DIMS variants."""
import pickle, itertools, collections, sys
exec(open('workshop/rounds/057/toolsmith_endt2.py').read())
from fractions import Fraction

def quotient_dims(arrows, rels):
    """dim of K Q/I per vertex pair, Q acyclic; paths = tuples of arrows in composition order (arrow = (s,h,k))."""
    out = collections.defaultdict(list)
    for a in arrows: out[a[0]].append(a)
    verts = sorted({a[0] for a in arrows} | {a[1] for a in arrows})
    paths = collections.defaultdict(list)       # (s,t) -> list of paths
    for v in verts: paths[(v, v)].append(())
    frontier = [((a,), a[0], a[1]) for a in arrows]
    L = 0
    while frontier:
        L += 1; assert L <= len(verts), 'cyclic quiver'
        nxt = []
        for p, s, t in frontier:
            paths[(s, t)].append(p)
            for a in out[t]: nxt.append((p + (a,), s, a[1]))
        frontier = nxt
    dims = {}
    for (s, t), ps in list(paths.items()):
        idx = {p: i for i, p in enumerate(ps)}; rows = []
        for r in rels:
            for (rp, ) in [(None,)]: pass
            r_s = next(iter(r))[0][0]; r_t = next(iter(r))[-1][1]
            for pre in paths.get((s, r_s), []):
                for post in paths.get((r_t, t), []):
                    row = {}
                    for rp, cf in r.items():
                        q = pre + tuple(rp) + post
                        if q in idx: row[idx[q]] = row.get(idx[q], 0) + Fraction(cf)
                    if row: rows.append(row)
        # rank
        piv = {}; rk = 0
        for row in rows:
            row = dict(row)
            while row:
                k = min(row)
                if k in piv:
                    f = row[k] / piv[k][k]
                    for kk, vv in piv[k].items():
                        row[kk] = row.get(kk, 0) - f * vv
                        if row[kk] == 0: del row[kk]
                else: piv[k] = row; rk += 1; break
        dims[(s, t)] = len(ps) - rk
    return dims

def lines_of(crels, P):
    """indices of relations of the form c0 (x,b0) + c1 (x,b1) with b0,b1 = P in the 2nd position and the same x (single or binomial)."""
    out = []
    for i, r in enumerate(crels):
        ks = list(r)
        if all(len(p) == 2 and p[1] in P for p in ks) and len({p[0] for p in ks}) == 1: out.append(i)
    return out

LINES = [(1, 0), (0, 1), (1, 1), (1, -1)]
def pattern(fn): return tuple(next(j for j in range(len(fn)) if fn[j] == fn[i]) for i in range(len(fn)))   # coincidence pattern (projective points)
def proj_eq(u, v): return u[0] * v[1] - u[1] * v[0] == 0

import hashlib; print('pickle sha256', hashlib.sha256(open(PK,'rb').read()).hexdigest()[:16], 'bytes', len(open(PK,'rb').read()))
verd = collections.Counter(); steps_tab = set(); steps_drop = set()
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
table = collections.Counter(); detail = []
for n_, r in enumerate(recs):
    a = r['parentObj']; v = r['v']; c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
    TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); A, d, carrows, crels = child_info(c)
    cc = collections.Counter((s, h) for (s, h, k) in carrows); par = [p for p, m in cc.items() if m > 1]
    if len(par) != 1 or cc[par[0]] != 2: continue
    P = sorted(ar for ar in carrows if (ar[0], ar[1]) == par[0]); L = lines_of(crels, P)
    if not L or len(L) > 3: continue
    # killed point of relation i: functional (c0,c1) on (x b0, x b1) = the line killed; as a point of P^1
    def pt(rl): return (rl.get((next(iter(rl))[0], P[0]), 0), rl.get((next(iter(rl))[0], P[1]), 0))
    true_pts = [pt(crels[i]) for i in L]
    xs = [next(iter(crels[i]))[0] for i in L]
    for choice in itertools.product(range(4), repeat=len(L)):
        new = list(crels)
        for i, ch, x in zip(L, choice, xs):
            lam = LINES[ch]; new[i] = {t: Fraction(w) for t, w in {(x, P[0]): lam[0], (x, P[1]): lam[1]}.items() if w != 0}
        pts = [LINES[ch] for ch in choice]
        qd = quotient_dims(carrows, new); td = {k: x for k, x in d.items() if x}
        eqdim = all(qd.get(k, 0) == d.get(k, 0) for k in set(qd) | set(d))
        truth = pattern(pts) == pattern(true_pts) if all(True for _ in pts) else None
        # pattern() uses exact tuple equality; make it projective
        def pat(ps): return tuple(next(j for j in range(len(ps)) if proj_eq(ps[j], ps[i])) for i in range(len(ps)))
        truth = (pat(pts) == pat(true_pts))
        vd = symcheck2(TE, arrows, 'same', carrows, new)
        verd[(('truth=NOT' if not truth else 'truth=iso'), vd)] += 1; steps_tab.add(n_)
        key = ('eqdim' if eqdim else 'DIMDIFF', 'truth=iso' if truth else 'truth=NOT', 'test=' + ('iso' if vd.startswith('iso') else 'NO'))
        table[key] += 1; detail.append((n_, v, len(L), choice, key))
    print('step', n_, 'v', v, 'k', len(L), 'true pts', true_pts, flush=True)
for k in sorted(table): print(k, table[k])
print('TOTAL variants', sum(table.values()))
print('VERDICT STRINGS (truth, verdict):'); [print('  ', k, c) for k, c in sorted(verd.items())]
print('steps in table', sorted(steps_tab))

# Part 2: dimension is an unchecked input of symcheck2.  Drop the first line-type relation of each step: dim K Q/I' > dim End(T), yet symcheck2 says iso.
dr = collections.Counter()
for n_, r in enumerate(recs):
    a = r['parentObj']; v = r['v']; c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
    TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); A, d, carrows, crels = child_info(c)
    cc = collections.Counter((s, h) for (s, h, k) in carrows)
    if not any(m > 1 for m in cc.values()): continue
    steps_drop.add(n_)
    for i in range(len(crels)):
        new = crels[:i] + crels[i + 1:]; qd = quotient_dims(carrows, new)
        big = sum(qd.values()) > sum(d.values()); vd = symcheck2(TE, arrows, 'same', carrows, new)
        dr[('dim(KQ/I\') > dim End(T)' if big else 'dim equal', 'test=' + ('iso' if vd.startswith('iso') else 'NO'))] += 1
print('DROP-ONE-RELATION (parallel-arrow steps):', dict(dr))

print('steps in drop-one', sorted(steps_drop))
