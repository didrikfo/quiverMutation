"""Lemma at l.1055-1059 of Oppermann 1504.02617v1: test d(dec(p) a^-1) = d(dec p) a^-1 + (-1)^|p| dec(p) a^* + s*(-1)^ci(t p) dec(dp/a)
for s=-1 (statement l.1058) and s=+1 (the l.1105 last line read literally as (-1)^{|p2|}.. see submission), on random paths p.
Also counts d^2 failures by letter type for the printed sign of the last term of d(a phi b^-1) (l.394).
Run: .venv/bin/python workshop/rounds/069/scholar_oppermann_lemma.py"""
import random, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from scholar_oppermann_sign import *

def paths_from_i(q, rng, maxlen=3):
    # random path (applied-first) starting at i not starting with an A arrow, length 1..maxlen
    for _ in range(50):
        v = q.i; w = []
        for _ in range(rng.randint(1, maxlen)):
            opts = [a for a in q.arr if q.arr[a][0] == v and not (a in q.A and not w)]
            if not opts: break
            a = rng.choice(opts); w.append(a); v = q.arr[a][1]
        if w and w[0] not in q.A and w[0] in q.arr: return tuple(w)
    return None

def test(sg, ntr=300):
    rng = random.Random(1); res = {}
    for _ in range(ntr):
        q = random_instance(rng); m = M(q, {'last': -1})
        w = paths_from_i(q, rng)
        if w is None: continue
        P = {w: 1}
        decP = m.dec(P)
        if not decP: continue
        tp = q.arr[w[-1]][1]; ci = 1 if tp == q.i else 0
        pdeg = sum(m.deg(L) for L in next(iter(decP)))   # ||p||
        for al in sorted(q.A):
            lhs = m.partial(m.attach_right(decP, al))
            rhs = m.attach_right(m.partial(decP), al)
            for u, c in {k: v for k, v in decP.items()}.items():
                add(rhs, (('s', al),) + u, c * (-1) ** pdeg)
            dp = q.dpath(w); T = {u[1:]: c for u, c in dp.items() if u[0] == al}
            for u, c in m.dec(T).items(): add(rhs, u, sg * (-1) ** ci * c)
            ok = (lhs == rhs)
            key = ('loop' if tp == q.i else 'nonloop-end', 'start-at-i-loop' if q.arr[w[0]][0] == q.arr[w[0]][1] else '')
            res.setdefault((key, ok), 0); res[(key, ok)] += 1
    return res

for sg in (-1, 1):
    print('coefficient s on (-1)^ci(t p) dec(dp/a):', sg, test(sg))
