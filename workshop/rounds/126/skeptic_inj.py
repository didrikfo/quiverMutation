"""Round 126 (skeptic), T11 / H-030: the 8 minimal length-3 parents where the gate over-admits (E-312).
(1) shape: for each failing (vertex, side) case, the source u where p |-> (p beta)_beta (F) / p |-> (beta p)_beta (R) is not injective,
    dim of e_u A e_i, the kernel, and whether it is EXACTLY E-080's shape: a square u->x->i, u->y->i, ONE arrow beta on the mutation side,
    kernel spanned by (u x i) - c (u y i), and the input relation (u x i beta) - c (u y i beta) (arrow length 3).
(2) injectivity test = Ladkani 1001.4765 Prop 2.3(c), two implementations: LIB = copy of mutation.isTiltingAtVertex on (quiver, relations)
    over Q (repo arithmetic, signed coefficients), OWN = mod 32003 on the round-118 engine; both against Hom-verified tilting
    (skeptic_convention.invariants: Hom(T,T[m]) = 0, m = +-1, +-2).
Pre-registered: HIT = injectivity says NO on all 8 failing cases and agrees with Hom on every case of the E-312 sample (N5 0-5999, N6 0-3999).
usage: skeptic_inj.py seeds | sweep KIND N SEED0 COUNT | planted        (from the repo root)"""
import sys, random
from fractions import Fraction
from collections import Counter
MYARGV = sys.argv[1:]
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/118'); sys.path.insert(0, 'workshop/rounds/091')
import skeptic_conv3 as c3
import skeptic_convention as sc
from quivermutation import arrowPaths as ap, procedure
P = sc.P
SEEDS = [('N', 5, 871)] + [('N', 6, s) for s in (514, 1629, 1806, 1944, 2257, 2794, 3277)]

def signed(c): return c - P if c > P // 2 else c
def repo_input_signed(A):                       # E-312 fix: -1 passed as -1, not 32002 (R-016)
    q, rels = c3.repo_input(A)
    return q, [{k: signed(v) for k, v in r.items()} for r in rels]

def injLib(quiver, relations, vertex):
    """copy of quivermutation.mutation.isTiltingAtVertex, taking (quiver, relations)"""
    outArrows = ap.arrowsOutOf(quiver, vertex)
    for i in quiver.nodes:
        if i == vertex: continue
        cands = ap.allPathsBetween(quiver, i, vertex)
        if not cands: continue
        qd = len(cands) - len(ap.idealBasis(quiver, relations, i, vertex))
        if qd == 0: continue
        rk, piv = 0, {}
        for path in cands:
            row = {}
            for beta in outArrows:
                red = ap.reduceAgainstPivots(ap.combination([path + (beta,)]), ap.idealBasis(quiver, relations, i, beta[1]))
                for k, v in red.items(): row[(beta, k)] = Fraction(v)
            while row:
                h = min(row)
                if h not in piv:
                    piv[h] = {k: v / row[h] for k, v in row.items()}; rk += 1; break
                f, pv = row[h], piv[h]
                row = {k: row.get(k, 0) - f * pv.get(k, 0) for k in set(row) | set(pv)}
                row = {k: v for k, v in row.items() if v != 0}
        if rk < qd: return False
    return True

def injOwn(A, i, side):
    """mod-P kernel of the one map, per source u; returns list of (u, basis, betas, kernel vectors)"""
    bad = []
    if side == 'F': betas = [k for k, (a, b, _) in enumerate(A.arrows) if a == i]
    else: betas = [k for k, (a, b, _) in enumerate(A.arrows) if b == i]
    for u in range(A.n):
        if u == i: continue
        basis = A.paths(u, i) if side == 'F' else A.paths(i, u)
        if not basis: continue
        cols = []; keys = {}
        for p in basis:
            col = {}
            for k in betas:
                img = A.compose(p, (k,)) if side == 'F' else A.compose((k,), p)
                for r, c in img.items():
                    key = keys.setdefault((k, r), len(keys)); col[key] = (col.get(key, 0) + c) % P
            cols.append(col)
        M = [[cols[j].get(r, 0) for j in range(len(basis))] for r in range(len(keys))]
        ker = sc.nullspace(M, len(basis)) if M else sc.nullspace([], len(basis))
        if ker: bad.append((u, basis, betas, ker))
    return bad

def vseq(A, p): return '-'.join(str(v + 1) for v in ([A.arrows[p[0]][0]] + [A.arrows[k][1] for k in p])) if p else '()'

def e080(A, i, side, bad):
    """classify a non-injective case. 'exact': one source u, one arrow beta on the side, dim e_u A e_i = 2, kernel spanned by
    p - c q (p, q nonzero length-2 paths u ~> i) for an INPUT relation (p - c q) beta (arrow length 3): E-080's shape.
    'core+k': the same kernel mechanism but dim e_u A e_i = 2 + k (extra paths u ~> i that survive). Otherwise 'other: <why>'."""
    if len(bad) != 1: return 'other: kernel at %d sources' % len(bad)
    u, basis, betas, ker = bad[0]
    if len(betas) != 1: return 'other: %d arrows on the side' % len(betas)
    if len(ker) != 1: return 'other: kernel dim %d' % len(ker)
    k = betas[0]; kv = ker[0]
    def coords(p):
        red = A.compose(p, ()) if p else {(): 1}
        return [red.get(b, 0) % P for b in basis]
    for r in A.rels:
        if len(r) != 2: continue
        (x, cx), (y, cy) = r.items()
        if len(x) != 3 or len(y) != 3: continue
        if side == 'F': ok = x[-1] == y[-1] == k; p, q = x[:-1], y[:-1]
        else: ok = x[0] == y[0] == k; p, q = x[1:], y[1:]
        if not ok: continue
        v = [(cx * a + cy * b) % P for a, b in zip(coords(p), coords(q))]
        if not any(v): continue
        if sc.rank([v, kv], len(basis)) == 1 and any(coords(p)) and any(coords(q)):
            tag = 'exact' if len(basis) == 2 else 'core+%d' % (len(basis) - 2)
            return '%s: %s / %s, arrow %s, coefs %s' % (tag, vseq(A, p), vseq(A, q), vseq(A, (k,)), [signed(cx), signed(cy)])
    return 'other: no input relation explains the kernel'

def homdegs(A, i, side):
    summ = [(v, sc.stalk(v)) for v in range(A.n) if v != i] + [(i, sc.mutate_cx(A, i, 'R' if side == 'F' else 'L'))]
    out = {}
    for m in (-2, -1, 1, 2):
        d = sum(sc.cx_hom(A, X, Y, m).dim for _, X in summ for _, Y in summ)
        if d: out[m] = d
    return out

def cases(A):
    q, rels = repo_input_signed(A)
    qo, ro = procedure._opposite(q, rels)
    for i in range(A.n):
        if any(a == b == i for a, b, _ in A.arrows): continue
        for side in ('F', 'R'):
            has = any((a if side == 'F' else b) == i for a, b, _ in A.arrows)
            if not has: continue
            gate = procedure.isMutable(q, rels, i + 1) if side == 'F' else procedure.isMutable(qo, ro, i + 1)
            lib = injLib(q, rels, i + 1) if side == 'F' else injLib(qo, ro, i + 1)
            own = injOwn(A, i, side)
            hom = sc.invariants(A, i, 'R' if side == 'F' else 'L')[0]
            yield i, side, gate, lib, own, hom

def seeds():
    for kind, n, seed in SEEDS:
        A = c3.gen(kind, n, random.Random(seed))
        print('== seed', seed, 'n', n, 'arrows', [(a + 1, b + 1) for a, b, _ in A.arrows],
              'rels', [{vseq(A, p): signed(c) for p, c in r.items()} for r in A.rels])
        for i, side, gate, lib, own, hom in cases(A):
            if gate and not hom:
                why = e080(A, i, side, own)
                print('  FAIL v%d %s gate %s Hom %s | LIB inj %s OWN inj %s | shape: %s | Hom(T,T[m]) %s' % (
                    i + 1, side, gate, hom, lib, not own, why, homdegs(A, i, side)))
                for u, basis, betas, ker in own:
                    print('     source %d: basis %s, arrows %s, kernel %s' % (u + 1, [vseq(A, p) for p in basis], [vseq(A, (k,)) for k in betas],
                          [{vseq(A, p): signed(c) for p, c in zip(basis, v) if c} for v in ker]))
            elif lib != hom or (not own) != hom:
                print('  DISAGREE v%d %s gate %s Hom %s LIB %s OWN %s' % (i + 1, side, gate, hom, lib, not own))

def sweep(kind, n, s0, cnt):
    tab = Counter(); dis = []
    for seed in range(s0, s0 + cnt):
        A = c3.gen(kind, n, random.Random(seed))
        if A is None: continue
        tab['parents'] += 1
        for i, side, gate, lib, own, hom in cases(A):
            tab['gate %s Hom %s LIB %s OWN %s' % (gate, hom, lib, not own)] += 1
            if lib != hom or (not own) != hom: dis.append((seed, i + 1, side, gate, hom, lib, not own))
            if gate and not hom:
                sh = e080(A, i, side, own); tab['gate over-admits; shape %s' % sh.split(':')[0]] += 1
                if sh.startswith('other'): print('  other: seed', seed, 'v%d' % (i + 1), side, sh)
    for k in sorted(tab): print(tab[k], k)
    print('injectivity != Hom:', len(dis), dis[:10])

def planted():
    """E-080's own algebra (n = 5, abde = acde, mutate at d, F): expect gate True, Hom False, inj False;
    abd = acd (length-2 square): expect inj True = Hom True; planted error: injLib with the relation sign flipped to +1 (abde + acde = 0
    is a different algebra but the same shape: still NO) and with the E-312 bug (-1 as 32002: x + 32002 y, a different algebra)."""
    a, b, c, d, e = range(5)
    arrows = [(a, b), (a, c), (b, d), (c, d), (d, e)]
    for name, rel in (('abde=acde', {(0, 2, 4): 1, (1, 3, 4): P - 1}), ('abd=acd', {(0, 2): 1, (1, 3): P - 1}),
                      ('abde+acde=0', {(0, 2, 4): 1, (1, 3, 4): 1})):
        A = c3.AlgG(5, arrows, [rel], True)
        for i, side, gate, lib, own, hom in cases(A):
            if i == d and side == 'F':
                print(name, 'v d F: gate %s Hom %s LIB %s OWN %s shape: %s' % (gate, hom, lib, not own, e080(A, i, side, own)))
    A = c3.AlgG(5, arrows, [{(0, 2, 4): 1, (1, 3, 4): P - 1}], True)
    q, rels = c3.repo_input(A)                  # unsigned: the R-016 harness bug, a different algebra
    print('abde=acde with -1 as 32002 (R-016 encoding): LIB at d', injLib(q, rels, d + 1), '(different algebra: x + 32002 y)')

if __name__ == '__main__':
    if MYARGV[0] == 'seeds': seeds()
    elif MYARGV[0] == 'planted': planted()
    else: sweep(MYARGV[1], int(MYARGV[2]), int(MYARGV[3]), int(MYARGV[4]))
