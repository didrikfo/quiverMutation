"""Round 107 (toolsmith), H-027 / T10: all indecomposable objects of K^b(proj B_n) up to shift, by CONE CLOSURE of the brute-force thin
complexes of toolsmith_enum.py (brute force alone finds 13/20 at n=5, 18/30 at n=6; D^b(B_n) ~ D^b(kD_n) has n(n-1) positive roots, checked
below by counting vectors v with v^T Cartan v = 1).  An object is kept iff End = k and Hom(X,X[m]) = 0 for m != 0 (exceptional => indecomposable).
Usage: toolsmith_cones.py n  -> workshop/rounds/107/toolsmith_objs_<n>.json  (list of [layers, d] with d entries [a,b,[path],coef])"""
import sys, itertools, json, time, random
sys.path.insert(0, '.')
exec(open('workshop/rounds/107/toolsmith_enum.py').read().split("if __name__")[0])
from fractions import Fraction

def rootsCount(n):
    from quivermutation import invariants
    C = np.array(invariants.cartanMatrix(Bn(n), exact=True).tolist(), dtype=int)
    return sum(1 for v in itertools.product(range(-4, 5), repeat=n) if np.array(v) @ C @ np.array(v) == 1) // 2

# complex := (degs {k:[vertices]}, d {k:{(a,b):{path:Fraction}}})
def shiftC(X, j):   # X[j]: degree k -> k - j  (no sign; isomorphic)
    degs, d = X
    return ({k - j: list(v) for k, v in degs.items()}, {k - j: {ab: dict(c) for ab, c in m.items()} for k, m in d.items()})

def normalize(X):
    degs, d = X; lo = min(k for k, v in degs.items() if v)
    return ({k - lo: v for k, v in degs.items() if v}, {k - lo: m for k, m in d.items() if m})

def chainMaps(A, X, Y):
    """basis of the space of chain maps X->Y (degree 0): list of {k:{(a,b):{path:coef}}}."""
    Xd, Xdi = X; Yd, Ydi = Y
    unk = []
    for k in sorted(set(Xd) & set(Yd)):
        for a, va in enumerate(Xd[k]):
            for b, vb in enumerate(Yd[k]):
                for p in A.basis(va, vb): unk.append((k, a, b, p))
    idx = {u: i for i, u in enumerate(unk)}
    rows = {}
    def add(key, u, c):
        rows.setdefault(key, {}); rows[key][idx[u]] = rows[key].get(idx[u], 0) + c
    for (k, a, b, p) in unk:
        va, vb = Xd[k][a], Yd[k][b]
        # d_Y f : X^k -> Y^{k+1}
        for (b1, c), cm in Ydi.get(k, {}).items():
            if b1 != b: continue
            vc = Yd[k + 1][c]
            comb = {p + pc: cc for pc, cc in cm.items()}
            for pp, cf in A.nf(va, vc, {q: Fraction(v) for q, v in comb.items()}).items(): add((k, a, c, pp), (k, a, b, p), cf)
        # f d_X : X^{k-1} -> Y^k   (f^k after d_X^{k-1})
        for (a0, a1), cm in Xdi.get(k - 1, {}).items():
            if a1 != a: continue
            v0 = Xd[k - 1][a0]
            comb = {pc + p: cc for pc, cc in cm.items()}
            for pp, cf in A.nf(v0, vb, {q: Fraction(v) for q, v in comb.items()}).items(): add((k - 1, a0, b, pp), (k, a, b, p), -cf)
    # nullspace of the system rows (sum_u rows[key][u] x_u = 0)
    M = [dict(r) for r in rows.values()]
    piv = {}
    for r in M:
        r = {u: Fraction(v) for u, v in r.items() if v}
        for h in sorted(piv):
            if h in r:
                f = r[h]
                for u, v in piv[h].items(): r[u] = r.get(u, 0) - f * v
                r = {u: v for u, v in r.items() if v}
        if r:
            h = min(r); inv = 1 / r[h]; r = {u: v * inv for u, v in r.items()}
            for h2 in list(piv):
                if h in piv[h2]:
                    f = piv[h2][h]
                    for u, v in r.items(): piv[h2][u] = piv[h2].get(u, 0) - f * v
                    piv[h2] = {u: v for u, v in piv[h2].items() if v}
            piv[h] = r
    free = [i for i in range(len(unk)) if i not in piv]
    basis = []
    for fcol in free:
        x = {fcol: Fraction(1)}
        for h, r in piv.items():
            if fcol in r: x[h] = -r[fcol]
        f = {}
        for i, v in x.items():
            k, a, b, p = unk[i]; f.setdefault(k, {}).setdefault((a, b), {})[p] = v
        basis.append(f)
    return basis

def cone(A, X, Y, f):
    Xd, Xdi = X; Yd, Ydi = Y
    ks = sorted(set(k - 1 for k in Xd) | set(Yd))
    layers = {}; off = {}
    for k in ks:
        layers[k] = list(Xd.get(k + 1, [])) + list(Yd.get(k, [])); off[k] = len(Xd.get(k + 1, []))
    d = {}
    for k in ks:
        m = {}
        if k + 1 in layers:
            for (a, b), c in Xdi.get(k + 1, {}).items(): m[(a, b)] = {p: -v for p, v in c.items()}          # -d_X : X^{k+1}->X^{k+2}
            for (a, b), c in f.get(k + 1, {}).items(): m[(a, off[k + 1] + b)] = dict(c)                      # f : X^{k+1} -> Y^{k+1}
            for (a, b), c in Ydi.get(k, {}).items(): m[(off[k] + a, off[k + 1] + b)] = dict(c)               # d_Y
        d[k] = m
    return (layers, d)

def minimize(A, C):
    layers, d = C
    layers = {k: list(v) for k, v in layers.items()}; d = {k: {ab: dict(c) for ab, c in m.items()} for k, m in d.items()}
    changed = True
    while changed:
        changed = False
        for k in sorted(d):
            for (a, b), c in list(d[k].items()):
                if layers[k][a] == layers[k + 1][b] and () in c and c[()] != 0:
                    inv = 1 / Fraction(c[()])
                    # eliminate a in layer k, b in layer k+1
                    new = {}
                    for (a2, b2), c2 in d[k].items():
                        if a2 == a or b2 == b: continue
                        new[(a2, b2)] = {p: v for p, v in c2.items()}
                    for (a2, bb), c1 in d[k].items():
                        if bb != b or a2 == a: continue
                        for (aa, b2), c3 in d[k].items():
                            if aa != a or b2 == b: continue
                            comb = {}
                            for p1, v1 in c1.items():
                                for p3, v3 in c3.items(): comb[p1 + p3] = comb.get(p1 + p3, 0) - v1 * inv * v3
                            comb = A.nf(layers[k][a2], layers[k + 1][b2], {q: Fraction(v) for q, v in comb.items()})
                            e = new.setdefault((a2, b2), {})
                            for p, v in comb.items(): e[p] = e.get(p, 0) + v
                    # reindex
                    def rm(lst, i): return lst[:i] + lst[i + 1:]
                    mapA = {i: (i if i < a else i - 1) for i in range(len(layers[k])) if i != a}
                    mapB = {i: (i if i < b else i - 1) for i in range(len(layers[k + 1])) if i != b}
                    d[k] = {(mapA[a2], mapB[b2]): {p: v for p, v in c.items() if v} for (a2, b2), c in new.items() if any(v for v in c.values())}
                    if k - 1 in d:   # maps into layer k: drop column a, reindex
                        d[k - 1] = {(x, mapA[y]): c for (x, y), c in d[k - 1].items() if y != a}
                    if k + 1 in d:   # maps out of layer k+1: drop row b
                        d[k + 1] = {(mapB[x], y): c for (x, y), c in d[k + 1].items() if x != b}
                    layers[k] = rm(layers[k], a); layers[k + 1] = rm(layers[k + 1], b)
                    changed = True; break
            if changed: break
    layers = {k: v for k, v in layers.items() if v}
    d = {k: {ab: {p: v for p, v in c.items() if v} for ab, c in m.items() if any(v for v in c.values())} for k, m in d.items() if k in layers and k + 1 in layers}
    return (layers, d)

def klassC(X, n):
    v = [0] * n
    for k, l in X[0].items():
        for x in l: v[x - 1] += (-1) ** k
    return tuple(v)
def keyOf(X, n):
    v = klassC(X, n); return max(v, tuple(-t for t in v))
def isExc(A, X):
    if homdim(A, X, X, 0) != 1: return False
    return not any(homdim(A, X, X, m) for m in (-4, -3, -2, -1, 1, 2, 3, 4))

def dsq(A, X):   # check d^2 = 0 (sanity)
    layers, d = X
    for k in d:
        if k + 1 not in d: continue
        res = {}
        for (a, b), c1 in d[k].items():
            for (b2, c), c2 in d[k + 1].items():
                if b2 != b: continue
                e = res.setdefault((a, c), {})
                for p1, v1 in c1.items():
                    for p2, v2 in c2.items(): e[p1 + p2] = e.get(p1 + p2, 0) + v1 * v2
        for (a, c), e in res.items():
            if A.nf(layers[k][a], layers[k + 2][c], {q: Fraction(v) for q, v in e.items()}): return False
    return True

def ser(X): return [sorted([k, v] for k, v in X[0].items()), sorted([k, sorted([list(ab), sorted([[list(p), str(v)] for p, v in c.items()])] for ab, c in m.items())] for k, m in X[1].items())]
def deser(S):
    degs = {k: v for k, v in S[0]}
    d = {k: {tuple(ab): {tuple(tuple(e) for e in p): Fraction(v) for p, v in c} for ab, c in m} for k, m in S[1]}
    return (degs, d)

if __name__ == '__main__':
    n = int(sys.argv[1]); random.seed(1)
    alg = Bn(n); A = Alg(alg); V = sorted(alg.quiver.nodes)
    print('n', n, 'positive roots of Euler form (= expected objects up to shift):', rootsCount(n)); sys.stdout.flush()
    t0 = time.time()
    cs = enumerate_complexes(A, V, 2, 5, 7)
    objs = {}
    for layers, ds in cs:
        X = normalize(cx(layers, ds))
        key = keyOf(X, n)
        if key in objs: continue
        if isExc(A, X): objs[key] = X
    print('brute force', len(objs), 'time %.0f' % (time.time() - t0)); sys.stdout.flush()
    target = rootsCount(n); rnd = 0
    while len(objs) < target and rnd < 6:
        rnd += 1; new = {}
        L = list(objs.values())
        for X in L:
            for Y in L:
                for j in range(-4, 5):
                    Yj = shiftC(Y, j)
                    if homdim(A, X, Yj, 0) == 0: continue
                    for f in chainMaps(A, X, Yj)[:3]:
                        Cm = minimize(A, cone(A, X, Yj, f))
                        if not Cm[0]: continue
                        Cn = normalize(Cm)
                        if not dsq(A, Cn): continue
                        key = keyOf(Cn, n)
                        if key in objs or key in new: continue
                        if isExc(A, Cn): new[key] = Cn
        objs.update(new)
        print('round', rnd, 'objects', len(objs), 'time %.0f' % (time.time() - t0)); sys.stdout.flush()
    json.dump([ser(X) for X in objs.values()], open('workshop/rounds/107/toolsmith_objs_%d.json' % n, 'w'))
    print('final', len(objs), 'of', target)
