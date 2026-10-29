"""Round 097 maverick: Z-congruence of the symmetrised Cartan form S = C + C^T (diag 2, even; C^-1 + C^-T = C^-1 S C^-T is Z-congruent to it)
for the certified-inequivalent profile-equal pairs of E-259.  Stage 1 (decided): signature, det, discriminant-form isomorphism q_S (brute force over the
finite group A = S^-1 Z^n / Z^n) => same genus (Nikulin 1.9.4 / 1.10.1); with rank >= l(A)+2 and indefinite, class = genus (Nikulin 1.13.3; from memory, UNVERIFIED).
Stage 2: Z[x]-module isomorphism (P Phi_a = Phi_b P, P in GL_n(Z)) searched by small combinations in the saturated commutant lattice (YES-only search).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/097/maverick_cartanZ.py N [members_per_quipu]"""
import sys, collections, itertools, random
from fractions import Fraction
import numpy as np
import sympy as sp
from quivermutation import coxeterTables as ct
from quivermutation.quipuForms import quipuForAlmostSeparateLNA
from cypari import pari

def cartan(N, rl):
    return sp.Matrix(np.array(ct.lnaCartanMatrix(N, rl)).astype(int).tolist())

def disc_form(S):
    n = S.shape[0]; d = int(S.det()); Si = S.inv()
    Sinv = [[Fraction(int(Si[i, j].p), int(Si[i, j].q)) for j in range(n)] for i in range(n)]
    def key(v): return tuple(sum(Sinv[i][j] * v[j] for j in range(n)) % 1 for i in range(n))
    def vec(k): # integer rep not needed; work in x-coordinates
        return k
    # elements as x = S^-1 v mod Z^n ; add generator columns
    gens = []; elems = {tuple([Fraction(0)] * n)}
    cols = [tuple(Sinv[i][j] % 1 for i in range(n)) for j in range(n)]
    def add(a, b): return tuple((x + y) % 1 for x, y in zip(a, b))
    for c in cols:
        if c in elems: continue
        gens.append(c); new = set(elems); cur = c
        while cur not in elems:
            new |= {add(e, cur) for e in elems}; cur = add(cur, c)
            if len(new) > 10**5: raise RuntimeError("group too big")
        elems = new
    Sf = [[Fraction(int(S[i, j])) for j in range(n)] for i in range(n)]
    def bil(x, y): return sum(x[i] * Sf[i][j] * y[j] for i in range(n) for j in range(n)) % 1
    def q(x): return (bil(x, x) / 2) % 1   # q(x)=x^T S x/2 mod 1 (S even)  (quadratic form values mod 1)
    return abs(d), elems, gens, bil, q, add

def order(x, add):
    z = tuple(Fraction(0) for _ in x); k = 1; c = x
    while c != z: c = add(c, x); k += 1
    return k

def disc_iso(Sa, Sb):
    da, Ea, Ga, bila, qa, adda = disc_form(Sa); db, Eb, Gb, bilb, qb, addb = disc_form(Sb)
    if da != db or len(Ea) != len(Eb): return False, "order"
    # q mod 2? S even: q takes values in Q/Z via x^T S x/2 mod 1 ; need Z/2Z-valued refinement: x^T S x mod 2 -> use qq = x^T S x mod 2
    n = Sa.shape[0]
    def qq(S, x):
        Sf = [[Fraction(int(S[i, j])) for j in range(n)] for i in range(n)]
        return sum(x[i] * Sf[i][j] * x[j] for i in range(n) for j in range(n))
    # x lives mod Z^n; x^T S x mod 2 is well defined only for x in S^-1 Z^n mod Z^n for even S: reduce reps mod 1 changes by even integer? yes (S even, x in dual).
    qa2 = lambda x: qq(Sa, x) % 2; qb2 = lambda x: qq(Sb, x) % 2
    # minimal generating set
    def gset(G, add):
        out = []; el = {tuple(Fraction(0) for _ in G[0])} if G else set()
        for g in G:
            if g in el: continue
            out.append(g); new = set(el); cur = g
            while cur not in el:
                new |= {add(e, cur) for e in el}; cur = add(cur, g)
            el = new
        return out
    ga = gset(Ga, adda); orders = [order(g, adda) for g in ga]
    candB = [[y for y in Eb if order(y, addb) == o and qb2(y) == qa2(g)] for g, o in zip(ga, orders)]
    for imgs in itertools.product(*candB):
        if any(bilb(imgs[i], imgs[j]) != bila(ga[i], ga[j]) for i in range(len(ga)) for j in range(i, len(ga))): continue
        mp = {}; ok = True
        for cs in itertools.product(*[range(o) for o in orders]):
            x = tuple(Fraction(0) for _ in range(n)); y = x
            for c, g, h in zip(cs, ga, imgs):
                for _ in range(c): x = adda(x, g); y = addb(y, h)
            if x in mp and mp[x] != y: ok = False; break
            mp[x] = y
        if ok and len(mp) == len(Ea) and len(set(mp.values())) == len(Eb) and all(qb2(mp[x]) == qa2(x) for x in mp): return True, imgs
    return False, "no isometry of finite quadratic forms"

def sig(S):
    ev = np.linalg.eigvalsh(np.array(S.tolist(), float)); return int((ev > 1e-9).sum()), int((ev < -1e-9).sum())

def phi(C): return -(C.T).inv() * C

def commutant_lattice(Pa, Pb):
    n = Pa.shape[0]
    # P Pa = Pb P  -> (Pa^T (x) I - I (x) Pb) vec(P) = 0 (row-major vec)
    A = sp.kronecker_product(Pa.T, sp.eye(n)) - sp.kronecker_product(sp.eye(n), Pb)
    M = pari.matrix(n * n, n * n, [int(v) for v in A])  # column-major in pari: build explicitly
    return pari.matkerint(M)

def module_iso_search(Pa, Pb, tries, seed=0):
    n = Pa.shape[0]
    A = np.kron(np.array(Pa.T.tolist(), dtype=object), np.eye(n, dtype=object)) - np.kron(np.eye(n, dtype=object), np.array(Pb.tolist(), dtype=object))
    m = n * n
    M = pari.matrix(m, m, [int(A[i, j]) for i in range(m) for j in range(m)])
    K = pari.matkerint(M); k = int(pari.matsize(K)[1])
    basis = [np.array([int(K[i, j]) for i in range(m)], dtype=object).reshape(n, n) for j in range(k)]
    rng = random.Random(seed); best = None
    for t in range(tries):
        cs = [rng.choice([-1, 0, 0, 1]) for _ in range(k)] if t else [1] + [0] * (k - 1)
        P = sum((c * B for c, B in zip(cs, basis)), np.zeros((n, n), dtype=object))
        d = int(pari.matdet(pari.matrix(n, n, [int(P[i, j]) for i in range(n) for j in range(n)])))
        if d != 0 and (best is None or abs(d) < best[0]): best = (abs(d), P)
        if abs(d) == 1: return True, k, 1, P
    return False, k, (best[0] if best else None), None

if __name__ == "__main__":
    N = int(sys.argv[1]); per = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
    for rl, k in keys.items(): groups[k].append(rl)
    for k, v in groups.items():
        q = collections.defaultdict(list)
        for rl in v:
            r = quipuForAlmostSeparateLNA(N, list(rl))
            if r is not None: q[r].append(rl)
        if len(q) < 2: continue
        names = list(q); reps = [(r, rl) for r in names for rl in q[r][:per]]
        data = []
        for r, rl in reps:
            C = cartan(N, rl); S = C + C.T; data.append((r, rl, C, S, int(S.det()), sig(S)))
        print("GROUP size", len(v), "quipus", names, "|det S|", sorted({abs(d[4]) for d in data}), "sig", sorted({d[5] for d in data}), flush=True)
        for (r1, a, C1, S1, d1, s1), (r2, b, C2, S2, d2, s2) in itertools.combinations(data, 2):
            if r1 == r2: tag = "within"
            else: tag = "ACROSS"
            ok, info = (disc_iso(S1, S2) if (d1 == d2 and s1 == s2) else (False, "det/sig"))
            print("  ", tag, a, b, "genus-equal(q_S iso)" if ok else "GENUS DIFFERS: " + str(info), flush=True)
