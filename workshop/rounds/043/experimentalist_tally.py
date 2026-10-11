"""Round 043 (experimentalist): n-tally for E-145 on the class walk (guarded as rounds/041/theorist_dump.py, or 'off' = key guard off to depth D as 042/theorist_dump_off.py).
For every gate-admitted step with J != 0 (J = perI): shape (|out v|, |supp J|, dim J, C'e_v = e_v, u = e_w - e_i), and on H1&H2 steps: s with F^s e_w = e_i (|s|<=8),
order of F, c_2 = (Y^T F^2)_ww, Q(x) = det(xC_B+C_B^T) - det(xC'+C'^T) lowest term (library C_B), and Q == x(adjS_ii-adjS_wi-adjS_iw).
Usage: experimentalist_tally.py n cls budget_sec [off DEPTH] [--plan]    (repo root)"""
import sys, time
from collections import Counter
import numpy as np, sympy as sp
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = [a for a in ARGV if a != '--plan']; plan = '--plan' in ARGV
n, cls, budget = int(_a[1]), int(_a[2]), float(_a[3])
off = len(_a) > 4 and _a[4] == 'off'; maxdepth = int(_a[5]) if off else 10**9
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
print('n', n, 'classes', [(i, len(classes[k])) for i, k in enumerate(order)][:6], 'chosen', cls, 'seeds', len(classes[base]))
if plan: sys.exit(0)
x = sp.symbols('x')
def pdet(C): return sp.Poly((x * sp.Matrix(C.tolist()) + sp.Matrix(C.T.tolist())).det(), x)
t0 = time.time(); seen = set(); frontier = []; nalg = 0; stop = False; depth = 0
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
cnt = Counter(); distinct = set(); distinct_all = set(); nrec = 0; tcache = {}
def tally(r):
    Jn = {i: j for i, j in r['J'].items() if j}
    if not Jn: return
    V, v, CA, CB = r['V'], r['v'], r['CA'], r['CB']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    Cp = rm @ CA @ rm.T; keep = [a for a in range(n) if a != k]
    cnt['J!=0 steps'] += 1
    qk = (CB.tobytes(), Cp.tobytes())
    if qk not in tcache:
        q = (pdet(CB) - pdet(Cp)).all_coeffs()[::-1]; low = next((j for j, a in enumerate(q) if a), None)
        tcache[qk] = (low, int(q[low]) if low is not None else None, sp.Poly(sum(c * x**j for j, c in enumerate(q)), x))
    low, lc, Qp = tcache[qk]
    cnt['Q lowest (deg, coeff)', (low, lc)] += 1
    colv = all(Cp[a, k] == 0 for a in keep); e = np.eye(n - 1, dtype=int)
    shape = (len(r['out']), len(Jn), tuple(sorted(Jn.values())), colv)
    if len(r['out']) != 1 or len(Jn) != 1 or not colv:
        cnt['OFF-SHAPE (H1 or |out|/|J| fails)', shape] += 1; distinct_all.add((CA.tobytes(), v)); return
    w = keep.index(V.index(r['out'][0])); i = keep.index(V.index(list(Jn)[0]))
    if not (Cp[k, keep] == e[w] - e[i]).all():
        cnt['OFF-SHAPE (u != e_w - e_i)', tuple(int(a) for a in Cp[k, keep])] += 1; distinct_all.add((CA.tobytes(), v)); return
    Z = Cp[np.ix_(keep, keep)]; key = (Z.tobytes(), w, i)
    cnt['H1&H2 steps'] += 1; distinct.add(key); distinct_all.add((CA.tobytes(), v))
    S = x * sp.Matrix(Z.tolist()) + sp.Matrix(Z.T.tolist()); A = S.adjugate()
    cnt['P1 formula == library Q', sp.expand(x * (A[i, i] - A[w, i] - A[i, w]) - Qp.as_expr()) == 0] += 1
    Y = np.rint(np.linalg.inv(Z)).astype(int); F = Z @ Y.T; Fi = np.rint(np.linalg.inv(F)).astype(int)
    P = np.eye(n - 1, dtype=int); ordF = None
    for t in range(1, 40):
        P = P @ F
        if (P == np.eye(n - 1, dtype=int)).all(): ordF = t; break
    def Pw(kk):
        M = np.eye(n - 1, dtype=int)
        for _ in range(abs(kk)): M = M @ (F if kk > 0 else Fi)
        return M
    ss = tuple(sg for sg in range(-8, 9) if (Pw(sg)[:, w] == e[i]).all())
    cnt['s set (|s|<=8)', ss] += 1; cnt['order of F', ordF] += 1
    cnt['c_2 = (Y^T F^2)_ww', int((Y.T @ Pw(2))[w, w])] += 1
    if 1 in ss: cnt['s=1: c_2', int(Y[:, w] @ F[:, i])] += 1
recs = []
while frontier and not stop and depth < maxdepth:
    depth += 1
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget: stop = True; break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        Q = alg.quiver; V = sorted(Q.nodes)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != V: continue
            nrec += 1
            J = perI(alg, v)
            if J:
                CA = np.array(invariants.cartanMatrix(alg), dtype=int); CB = np.array(invariants.cartanMatrix(ch), dtype=int)
                tally(dict(CA=CA, CB=CB, V=V, v=v, J=dict(J), out=[w for (a, w) in Q.edges() if a == v]))
            if off or search._coxeterKeyOrNone(ch) == base: add(ch)
print('n', n, 'cls', cls, 'mode', 'off%d' % maxdepth if off else 'guarded', 'expanded', nalg, 'gate-admitted records', nrec, 'depth', depth, 'stopped(cap hit)', stop, 'secs %.0f' % (time.time() - t0))
print('distinct (Z,w,i) on H1&H2:', len(distinct), ' distinct (C_A,v) over all J!=0:', len(distinct_all))
for kx in sorted(cnt, key=str): print(kx, cnt[kx])
