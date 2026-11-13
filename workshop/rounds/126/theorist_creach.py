"""Round 126 (theorist), T5 / H-023: is condition (C) of E-322 reachable?  A parent where (C) -- and (LF) -- fail for EVERY
positive homogeneous weight of one presentation, and hold after a change of arrows (same algebra, isomorphic model).
Q: a1: i->1, c, e: 1->2 (parallel), b: i->2, x: 2->3, y: 2->4   (vertices i = 0, 1, 2, 3, 4; acyclic)
Presentation P:  R1 = b x - a1 c x,  R2 = b y - a1 e y,  R3 = a1 c + a1 e.
  Homogeneity forces w(b) = w(a1) + w(c) = w(a1) + w(e) = w(R3), so (C) (w(R3) > w(b)) fails for every weight, and the kM letter
  R3/b : 2 -> 2 (degree 1) has weight w(R3) - w(b) = 0: a weight-0 cycle, (LF) NO for every weight.
Presentation P': f = c + e, b' = b + a1 e  ->  R1' = b' x,  R2' = b' y - 2 a1 e y,  R3' = a1 f.  Here w(f) is free; with
  w(f) = 2 (C) holds and the (C) shift makes every letter positive.
Checks: [ISO] the substitution c = f - e, b = b' - a1 e maps ideal(P) onto ideal(P') (path spaces, exact); [MM] both models have
H = A cell by cell (exact rationals; no letters of degree >= 2 are needed); [LF] Floyd-Warshall on the kM letter quiver at i
for several weight vectors; [C-shift] min shifted letter weight.
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/126/theorist_creach.py"""
import sys, importlib.util, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '.')
def load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
P = load('lnap', 'workshop/rounds/123/theorist_lnaparent.py'); h = P.h; lf = P.lf; opp = P.opp
ARR0 = {'a1': (0, 1, 0), 'c': (1, 2, 0), 'e': (1, 2, 0), 'b': (0, 2, 0), 'x': (2, 3, 0), 'y': (2, 4, 0)}
REL_P = {'R1': {('b', 'x'): 1, ('a1', 'c', 'x'): -1}, 'R2': {('b', 'y'): 1, ('a1', 'e', 'y'): -1}, 'R3': {('a1', 'c'): 1, ('a1', 'e'): 1}}
ARR1 = {'a1': (0, 1, 0), 'f': (1, 2, 0), 'e': (1, 2, 0), 'b2': (0, 2, 0), 'x': (2, 3, 0), 'y': (2, 4, 0)}
REL_P1 = {'R1': {('b2', 'x'): 1}, 'R2': {('b2', 'y'): 1, ('a1', 'e', 'y'): -2}, 'R3': {('a1', 'f'): 1}}
SRC = {'R1': (0, 3), 'R2': (0, 4), 'R3': (0, 2)}

def paths(arr, s, t):
    out = []
    def rec(v, w):
        if v == t and w: out.append(tuple(w))
        for L, (a, b_, _) in arr.items():
            if a == v: rec(b_, w + [L])
    rec(s, []); return out

def ideal_span(arr, rels, s, t):
    """span of u r v (u: s -> s(r), v: t(r) -> t) inside the path space s -> t"""
    B = {}
    for R, E in rels.items():
        rs, rt = SRC[R]
        for u in ([()] if rs == s else [p for p in paths(arr, s, rs)]):
            for v in ([()] if rt == t else [p for p in paths(arr, rt, t)]):
                P.rref_add(B, {u + w + v: c for w, c in E.items()})
    return B

def dimA(arr, rels, s, t):
    return len(paths(arr, s, t)) - len(ideal_span(arr, rels, s, t)) + (1 if s == t else 0)

def subst(E):   # P -> P' generators: c = f - e, b = b2 - a1 e
    rule = {'c': {('f',): 1, ('e',): -1}, 'b': {('b2',): 1, ('a1', 'e'): -1}}
    out = {}
    for w, k in E.items():
        cur = {(): Fr(k)}
        for x in w:
            img = rule.get(x, {(x,): 1}); new = {}
            for u, c in cur.items():
                for v, e in img.items(): new[u + v] = new.get(u + v, 0) + c * e
            cur = {u: c for u, c in new.items() if c}
        for u, c in cur.items(): out[u] = out.get(u, 0) + c
    return {u: c for u, c in out.items() if c}

def model(arr, rels):
    a = dict(arr); d = {}
    for R, E in rels.items(): a[R] = (*SRC[R], 1); d[R] = dict(E)
    return a, d

def check_mm(arr, rels, model_rels=None):
    a, d = model(arr, model_rels or rels); q = opp.Q(a, [], 0); q.d = d; bad = []
    for s in range(5):
        for t in range(5):
            if s == t: continue
            ws = paths(a, s, t); bydeg = {}
            for w in ws: bydeg.setdefault(q.deg(w), []).append(w)
            for m in range(0, 3):
                z = len(P.kernel(bydeg.get(m, []), [q.dpath(w) for w in bydeg.get(m, [])])) if m else len(bydeg.get(0, []))
                B = {}
                for w in bydeg.get(m + 1, []): P.rref_add(B, q.dpath(w))
                want = dimA(arr, rels, s, t) if m == 0 else 0
                if z - len(B) != want: bad.append((s, t, m, z - len(B), want))
    return bad

def lf_at_i(arr, rels, w):
    a, d = model(arr, rels)
    A = sorted(L for L, (s, t, g) in a.items() if s == 0 and g == 0)
    q = lf.mk(a, A, d, w); q.i = 0
    m = opp.M(q, {'last': -1}); letters = opp.all_letters(q) + [('s', al) for al in A]
    wt = lambda L: -q.w[L[1]] if L[0] == 's' else q.w[L[2]] + (q.w[L[1]] if L[1] else 0) - (q.w[L[3]] if L[3] else 0)
    d2 = sum(1 for L in letters if m.partial(m.partial({(L,): 1})))
    ok, mins = P.lf_line(m, letters, wt)
    zero = [(h.nm(L), m.deg(L)) for L in letters if L[0] == 'm' and wt(L) <= 0]
    mA = max(w[x] for x in A); mN = min(w[R] for R in rels)
    C = mA < mN; sh = None
    if C:
        c = Fr(mA + mN, 2); phi = lambda v: c if v == 0 else 0
        sh = min(wt(L) + phi(m.tgt(L)) - phi(m.src(L)) for L in letters)
    return dict(d2=d2, LF=ok, C=C, nonpositive_letters=zero, min_shift=sh)

def weights(arr, rels, free):
    w = dict(free)
    for R, E in rels.items(): w[R] = sum(w[x] for x in next(iter(E)))
    homog = all(len({sum(w[x] for x in u) for u in E}) == 1 for E in rels.values())
    return w, homog

if __name__ == '__main__':
    # [ISO] substitution maps ideal(P) onto ideal(P'), cell by cell
    iso = []
    for s, t in itertools.product(range(5), repeat=2):
        if s == t: continue
        I0 = ideal_span(ARR0, REL_P, s, t); I1 = ideal_span(ARR1, REL_P1, s, t)
        img = {}
        for row in I0.values(): P.rref_add(img, subst(row))
        both = dict(I1)
        for row in img.values(): P.rref_add(both, row)
        if not (len(img) == len(I1) == len(both)): iso.append((s, t, len(img), len(I1), len(both)))
    print(f'[ISO] c = f - e, b = b2 - a1 e maps ideal(P) onto ideal(P\'): {"YES" if not iso else "NO " + str(iso)}')
    iso_ctrl = subst({('b', 'x'): 1})   # control: R1 without its a1 c x term is not in ideal(P')
    I1 = ideal_span(ARR1, REL_P1, 0, 3); t_ = dict(I1); grows = P.rref_add(t_, iso_ctrl)
    print(f'[ISO control] image of b x alone lies in ideal(P\'): {not grows} (must be False)')
    print('dim e_s A e_t:', {(s, t): dimA(ARR0, REL_P, s, t) for s in range(5) for t in range(5) if s != t and dimA(ARR0, REL_P, s, t)})
    for name, arr, rels in (('P', ARR0, REL_P), ("P'", ARR1, REL_P1)):
        print(f'[MM] {name}: H = A in every cell (degrees 0-2): {"YES" if not check_mm(arr, rels) else "NO " + str(check_mm(arr, rels)[:3])}')
    bad = {k: v for k, v in REL_P.items() if k != 'R2'}   # control: relation letter R2 left out of the model
    print(f'[MM control] P model without the letter R2: failures {len(check_mm(ARR0, REL_P, bad))} (must be > 0; a dropped TERM of d R1 is'
          f' NOT caught -- {len(check_mm(ARR0, REL_P, dict(REL_P, R1={("b", "x"): 1})))} failures -- the check compares dimensions only)')
    for free in ({'a1': 1, 'c': 1, 'e': 1, 'b': 2, 'x': 1, 'y': 1}, {'a1': 2, 'c': 3, 'e': 3, 'b': 5, 'x': 1, 'y': 4}, {'a1': 1, 'c': 5, 'e': 5, 'b': 6, 'x': 2, 'y': 1}):
        w, ok = weights(ARR0, REL_P, free)
        print(f'[P]  weights {free} homogeneous {ok}: {lf_at_i(ARR0, REL_P, w)}')
    w, ok = weights(ARR0, REL_P, {'a1': 1, 'c': 1, 'e': 2, 'b': 2, 'x': 1, 'y': 1})
    print(f'[P]  non-homogeneous control w(e) = 2: homogeneous {ok} (must be False: w(b) = w(a1) + w(c) = w(a1) + w(e) is forced)')
    for free in ({'a1': 1, 'f': 2, 'e': 1, 'b2': 2, 'x': 1, 'y': 1}, {'a1': 1, 'f': 1, 'e': 1, 'b2': 2, 'x': 1, 'y': 1}):
        w, ok = weights(ARR1, REL_P1, free)
        print(f"[P'] weights {free} homogeneous {ok}: {lf_at_i(ARR1, REL_P1, w)}")
