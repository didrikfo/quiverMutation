"""Round 107 (toolsmith), H-027 / T10: enumerate ALL indecomposable objects of K^b(proj B_n) up to shift (exceptional minimal complexes
with thin layers), using the round-050 chain-map Hom (skeptic_tilt.homdim).  Completeness check: D^b(A_n)-type has n(n+1)/2 orbits.
Usage: toolsmith_enum.py [n] [maxTermsPerLayer] [maxLen]   -> writes /tmp-free pickle path given by env ENUMOUT (default workshop/rounds/107/enum_n.json)"""
import sys, itertools, json, os, time
sys.path.insert(0, '.')
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import numpy as np
from quivermutation import nakayama as nk

def LNA(n, row): return nk.LinearNakayamaAlgebra(n, list(row))
def Bn(n): return LNA(n, (0,) * (n - 4) + (3, 0))
def An(n): return LNA(n, (n - 1,) + (0,) * (n - 3))

def cx(layers, ds, S=0):
    """layers: list of tuples of vertices (degree S+i), ds: list of dict {(a,b): (path,coef)}"""
    degs = {S + i: list(l) for i, l in enumerate(layers)}
    d = {S + i: {ab: {pc[0]: Fraction(pc[1])} for ab, pc in ds[i].items()} for i in range(len(ds))}
    return (degs, d)

def d2zero(A, layers, dprev, dnext):
    # (d_next o d_prev): layer i -> i+2 ; entries compose paths; reduce mod I
    L0, L1, L2 = layers
    res = {}
    for (a, b), (p1, c1) in dprev.items():
        for (b2, c), (p2, c2) in dnext.items():
            if b2 != b: continue
            res.setdefault((a, c), {})
            pth = p2 + p1   # same convention as compose_mat: g after f = pg + pf
            res[(a, c)][pth] = res[(a, c)].get(pth, 0) + c1 * c2
    for (a, c), comb in res.items():
        if A.nf(L0[a], L2[c], {k: Fraction(v) for k, v in comb.items()}): return False
    return True

def enumerate_complexes(A, V, maxPer, maxLen, maxTerms):
    layersets = [t for r in range(1, maxPer + 1) for t in itertools.combinations(V, r)]
    def allowed(La, Lb):
        out = []
        for a, x in enumerate(La):
            for b, y in enumerate(Lb):
                if x != y:
                    bs = A.basis(x, y)
                    if bs: out.append(((a, b), bs[0]))
        return out
    def layerMaps(La, Lb):
        al = allowed(La, Lb)
        for sel in itertools.product((0, 1, -1), repeat=len(al)):
            yield {ab: (p, c) for (ab, p), c in zip(al, sel) if c}
    res = []
    def rec(layers, ds, terms):
        if len(layers) >= 1 and (len(layers) == 1 and len(layers[0]) == 1 or ds): res.append((list(layers), list(ds)))
        if len(layers) >= maxLen: return
        for Lb in layersets:
            if terms + len(Lb) > maxTerms: continue
            for dm in layerMaps(layers[-1], Lb):
                if not dm: continue
                if len(layers) >= 2 and not d2zero(A, layers[-2:] + [Lb], ds[-1], dm): continue
                layers.append(Lb); ds.append(dm)
                rec(layers, ds, terms + len(Lb))
                layers.pop(); ds.pop()
    for L0 in layersets:
        rec([L0], [], len(L0))
    return res

def klass(layers, nV):
    v = [0] * nV
    for i, l in enumerate(layers):
        for x in l: v[x - 1] += (-1) ** i
    return tuple(v)

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    maxPer = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    maxLen = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    maxTerms = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    alg = Bn(n); A = Alg(alg); V = sorted(alg.quiver.nodes)
    t0 = time.time()
    cs = enumerate_complexes(A, V, maxPer, maxLen, maxTerms)
    print('candidates', len(cs), 'time %.1f' % (time.time() - t0)); sys.stdout.flush()
    good = {}
    for layers, ds in cs:
        X = cx(layers, ds)
        e = homdim(A, X, X, 0)
        if e != 1: continue
        if any(homdim(A, X, X, m) for m in (-3, -2, -1, 1, 2, 3)): continue
        key = klass(layers, n)
        key = max(key, tuple(-t for t in key))
        if key not in good: good[key] = (layers, {str(k): [v[0], v[1]] for k, v in enumerate(ds)} if False else [[(list(ab), [list(p[0]), p[1]]) for ab, p in d.items()] for d in ds])
    print('exceptional classes (up to sign)', len(good), 'expected', n * (n + 1) // 2, 'time %.1f' % (time.time() - t0))
    out = os.environ.get('ENUMOUT', 'workshop/rounds/107/toolsmith_enum_%d.json' % n)
    json.dump([[k, v[0], v[1]] for k, v in good.items()], open(out, 'w'))
