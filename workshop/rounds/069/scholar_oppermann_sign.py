"""d^2 test of Oppermann arXiv:1504.02617v1 Def 3.2 (l. 383-395) and Lemma at l. 1055-1059.
Words are tuples of letters in APPLIED-FIRST order. Leibniz: d(x_m..x_1)=sum_k (-1)^{sum_{l>k}deg x_l} ...
Usage: .venv/bin/python workshop/rounds/069/scholar_oppermann_sign.py [ntrials]
"""
import random, sys, itertools
from collections import defaultdict

def add(D, w, c):
    D[w] = D.get(w, 0) + c
    if D[w] == 0: del D[w]

class Q:
    def __init__(s, arrows, A_names, i):
        s.arr = arrows            # name -> (src, tgt, deg)
        s.A = set(A_names); s.i = i; s.d = {}
    def deg(s, w): return sum(s.arr[x][2] for x in w)
    def dpath(s, w):
        out = {}
        for k, x in enumerate(w):
            sign = (-1) ** sum(s.arr[y][2] for y in w[k+1:])
            for u, c in s.d.get(x, {}).items():
                add(out, w[:k] + u + w[k+1:], sign * c)
        return out
    def dlin(s, E):
        out = {}
        for w, c in E.items():
            for u, e in s.dpath(w).items(): add(out, u, c * e)
        return out

# M letters: ('m', a, phi, b) = a phi b^{-1} (a,b None or A-name); ('s', alpha) = alpha^*
class M:
    def __init__(s, q, sg_loop=None):
        s.q = q; s.i = q.i; s.sg = sg_loop or {}
    def src(s, L):
        q = s.q
        if L[0] == 's': return q.arr[L[1]][1]
        _, a, p, b = L
        return q.arr[b][1] if b else q.arr[p][0]
    def tgt(s, L):
        q = s.q
        if L[0] == 's': return s.i
        _, a, p, b = L
        return q.arr[a][1] if a else q.arr[p][1]
    def deg(s, L):
        q = s.q
        if L[0] == 's': return 0
        _, a, p, b = L
        sp, tp, dp = q.arr[p]
        d = dp - (sp == s.i) + (tp == s.i)
        return d - (1 if a else 0) + (1 if b else 0)
    def wdeg(s, w): return sum(s.deg(L) for L in w)
    # --- path conversion: Q-path (applied first) -> M element, with Delta expansion
    def dec_path(s, w, red=True):
        q = s.q
        if not w: return {(): 1}
        if red and w[0] in q.A: return {}
        letters = []
        k = 0
        while k < len(w):
            x = w[k]
            if x in q.A:   # merge into previous letter (previous is plain phi)
                prev = letters.pop()
                assert prev[1] is None and prev[3] is None
                letters.append(('m', x, prev[2], None))
            else: letters.append(('m', None, x, None))
            k += 1
        # Delta at junctions
        res = {tuple(letters): 1}
        for jx in range(len(letters) - 1):
            pass
        out = {}
        junc = [jx for jx in range(len(letters) - 1)
                if s.tgt(letters[jx]) == s.i and letters[jx][1] is None
                and letters[jx+1][3] is None and s.src(letters[jx+1]) == s.i
                and letters[jx+1][2] not in q.A]
        for choice in itertools.product([None] + sorted(q.A), repeat=len(junc)):
            Ls = list(letters); c = 1
            for jx, be in zip(junc, choice):
                if be is None: continue
                c *= -1
                a = Ls[jx]; b = Ls[jx+1]
                Ls[jx] = ('m', be, a[2], a[3]); Ls[jx+1] = ('m', b[1], b[2], be)
            add(out, tuple(Ls), c)
        return out
    def dec(s, E):
        out = {}
        for w, c in E.items():
            for u, e in s.dec_path(w).items(): add(out, u, c * e)
        return out
    def attach_right(s, E, be):   # E * be^{-1}: attach to first-applied letter
        out = {}
        for w, c in E.items():
            L = w[0]; assert L[0] == 'm' and L[3] is None
            add(out, (('m', L[1], L[2], be),) + w[1:], c)
        return out
    def partial_letter(s, L):
        q = s.q; i = s.i
        if L[0] == 's': return {}
        _, a, p, b = L
        sp, tp, dp = q.arr[p]
        out = {}
        def red_dphi():
            return {w: c for w, c in q.d.get(p, {}).items() if w[0] not in q.A}
        if a is None and b is None:
            for w, c in s.dec(red_dphi()).items(): add(out, w, -c if tp == i else c)
            if tp == i:
                for al in q.A: add(out, (('m', al, p, None), ('s', al)), 1)
            return out
        if a is not None and b is None:        # alpha phi
            E = {}
            for w, c in red_dphi().items(): add(E, w + (a,), c)
            return s.dec(E) if True else {}
        # b not None
        if a is None:
            base = s.partial_letter(('m', None, p, None))
            for w, c in s.attach_right(base, b).items(): add(out, w, c)
            add(out, (('s', b), ('m', None, p, None)), (-1) ** s.deg(('m', None, p, None)))
            sg = -((-1) ** (tp == i))
            T = {w[1:]: c for w, c in q.d.get(p, {}).items() if w[0] == b}
            for w, c in s.dec(T).items(): add(out, w, sg * c)
            return out
        base = s.partial_letter(('m', a, p, None))
        for w, c in s.attach_right(base, b).items(): add(out, w, c)
        add(out, (('s', b), ('m', a, p, None)), (-1) ** s.deg(('m', a, p, None)))
        T = {w[1:]: c for w, c in q.d.get(p, {}).items() if w[0] == b}
        Ta = {}
        for w, c in T.items(): add(Ta, w + (a,), c)
        for w, c in s.dec(Ta).items(): add(out, w, s.sg.get('last', 1) * c)
        return out
    def partial(s, E):
        out = {}
        for w, c in E.items():
            for k, L in enumerate(w):
                sign = (-1) ** s.wdeg(w[k+1:])
                for u, e in s.partial_letter(L).items():
                    add(out, w[:k] + u + w[k+1:], c * sign * e)
        return out

def random_instance(rng, deg2=True):
    # vertices 0=i,1,2,3
    arrows = {}
    n = 0
    def new(sv, tv, dg):
        nonlocal n; n += 1; arrows['x%d' % n] = (sv, tv, dg); return 'x%d' % n
    A = [new(0, 1, 0)]
    if rng.random() < .5: A.append(new(0, 2, 0))
    if rng.random() < .5: A.append(new(0, 1, 0))
    for _ in range(rng.randint(1, 2)): new(rng.choice([1, 2, 3]), 0, 0)   # into i, deg 0
    new(1, 3, 0); new(3, 2, 0)
    for _ in range(rng.randint(1, 3)): new(0, 0, 1)                       # loops deg1
    for _ in range(rng.randint(1, 3)): new(rng.randint(0, 3), rng.randint(0, 3), 1)
    q = Q(arrows, A, 0)
    deg0 = [a for a in arrows if arrows[a][2] == 0]
    def paths0(sv, tv, maxlen=4):
        out = []
        def rec(v, w):
            if v == tv and w: out.append(tuple(w))
            if len(w) < maxlen:
                for a in deg0:
                    if arrows[a][0] == v: rec(arrows[a][1], w + [a])
        rec(sv, []); return out
    deg1 = [a for a in arrows if arrows[a][2] == 1]
    for a in deg1:
        sv, tv, _ = arrows[a]
        P = paths0(sv, tv)
        E = {}
        for w in rng.sample(P, min(len(P), rng.randint(0, 3))): add(E, w, rng.choice([1, -1, 2]))
        q.d[a] = E
    if deg2:
        for _ in range(2):
            a, b = rng.choice(deg1), rng.choice(deg1)
            if arrows[a][0] == arrows[b][1]:   # b first then a
                n += 1; nm = 'x%d' % n
                arrows[nm] = (arrows[b][0], arrows[a][1], 2)
                E = {}
                for u, c in q.d[a].items(): add(E, (b,) + u, c)       # da * b
                for u, c in q.d[b].items(): add(E, u + (a,), -c)      # a * db
                q.d[nm] = E
    return q

def all_letters(q):
    L = []
    for p, (sv, tv, dg) in q.arr.items():
        if p in q.A: continue
        for a in [None] + (sorted(q.A) if tv == q.i else []):
            for b in [None] + (sorted(q.A) if sv == q.i else []):
                L.append(('m', a, p, b))
    return L

def run(sg, ntr, seed=0):
    rng = random.Random(seed); bad = 0; tot = 0; dbad = 0
    for _ in range(ntr):
        q = random_instance(rng)
        for a in q.arr:                      # sanity d^2 = 0 on kQ
            if q.dlin(q.d.get(a, {})): dbad += 1
        m = M(q, sg)
        for L in all_letters(q):
            tot += 1
            if m.partial(m.partial({(L,): 1})): bad += 1
    return dbad, bad, tot

if __name__ == '__main__':
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    for sgn in (1, -1):
        print('sign of last term of d(a phi b^-1) =', sgn, ' (d^2_kQ failures, d^2_M failures, letters):', run({'last': sgn}, ntr))
