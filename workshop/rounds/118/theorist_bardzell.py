"""Round 118 (theorist): independent check of HH^0..2 of quadratic monomial algebras by the Bardzell/Koszul complex (NOT the bar complex of
rounds/115/theorist_hh.py).  A = kQ/(paths a b in Rel), Q acyclic.  Minimal bimodule resolution: P_n = sum over AP_n = {a1..an : each a_i a_{i+1} in Rel}
(AP_0 = vertices, AP_1 = arrows) of A e_t (x) e_s A.  Cochains C^n = {(p, q): p in AP_n, q a nonzero path of A parallel to p};
(d phi)(a1..a_{n+1}) = a1 phi(a2..a_{n+1}) + (-1)^{n+1} phi(a1..an) a_{n+1}   (n >= 1);  n = 0: (d phi)(a) = a phi(e_t(a)) - phi(e_s(a)) a.
Path composition: a b = 'a then b'.  Exact integer rank over Q (sympy)."""
import sympy, itertools
def hh(V, arrows, rel, top=3):
    # arrows: name -> (s, t); rel: set of (a, b) with a then b = 0.  A-basis: nonzero paths (tuples of arrow names), vertices as ('e', v)
    S = {a: arrows[a][0] for a in arrows}; T = {a: arrows[a][1] for a in arrows}
    def nz(p): return all((p[i], p[i + 1]) not in rel for i in range(len(p) - 1))
    paths = [(('e', v),) for v in V]; cur = [(a,) for a in arrows]
    while cur:
        paths += cur; cur = [p + (a,) for p in cur for a in arrows if T[p[-1]] == S[a] and nz(p + (a,))]
    src = lambda p: p[0][1] if p[0][0] == 'e' else S[p[0]]
    tgt = lambda p: p[0][1] if p[0][0] == 'e' else T[p[-1]]
    def mul(p, q):   # p then q, None if zero
        if tgt(p) != src(q): return None
        if p[0][0] == 'e': return q
        if q[0][0] == 'e': return p
        r = p + q; return r if nz(r) else None
    AP = {0: [(('e', v),) for v in V], 1: [(a,) for a in arrows]}
    for n in range(2, top + 2):
        AP[n] = [c + (a,) for c in AP[n - 1] for a in arrows if T[c[-1]] == S[a] and (c[-1], a) in rel] if n > 1 else []
    ends = lambda c: (src(c), tgt(c))
    CB = {n: [(c, q) for c in AP[n] for q in paths if (src(q), tgt(q)) == ends(c)] for n in AP}
    ix = {n: {x: i for i, x in enumerate(CB[n])} for n in CB}
    def D(n):
        M = sympy.zeros(len(CB[n]), len(CB[n + 1]))
        for j, (c, q) in enumerate(CB[n + 1]):
            pass
        for (ch, tq) in CB[n + 1]:
            pass
        for ch in AP[n + 1]:
            if n == 0:
                a = ch[0]
                for (c0, b), i in ix[0].items():
                    v = c0[0][1]
                    if v == T[a]:
                        r = mul((a,), b)
                        if r is not None: M[i, ix[1][(ch, r)]] += 1
                    if v == S[a]:
                        r = mul(b, (a,))
                        if r is not None: M[i, ix[1][(ch, r)]] -= 1
                continue
            for (c2, b), i in ix[n].items():
                if c2 == ch[1:]:
                    r = mul((ch[0],), b)
                    if r is not None: M[i, ix[n + 1][(ch, r)]] += 1
                if c2 == ch[:-1]:
                    r = mul(b, (ch[-1],))
                    if r is not None: M[i, ix[n + 1][(ch, r)]] += (-1) ** (n + 1)
        return M
    rk = {n: D(n).rank() if len(CB[n]) and len(CB[n + 1]) else 0 for n in range(top + 1)}
    return [len(CB[n]) - rk[n] - (rk[n - 1] if n else 0) for n in range(top + 1)], {n: len(CB[n]) for n in CB}
if __name__ == '__main__':
    A3 = {'a': (1, 2), 'b': (2, 3)}
    print('A3 path algebra           ', hh([1, 2, 3], A3, set()), ' expect [1,0,0,..]')
    print('A3 / rad^2                ', hh([1, 2, 3], A3, {('a', 'b')}), ' expect [1,0,0,..]')
    print('Kronecker                 ', hh([1, 2], {'a': (1, 2), 'b': (1, 2)}, set()), ' expect [1,3,0,..]')
    sq = {'a': (1, 2), 'b': (2, 4), 'c': (1, 3), 'd': (3, 4)}
    print('square, ab = 0            ', hh([1, 2, 3, 4], sq, {('a', 'b')}), ' bar-complex code (round 115): (1,1,1)')
    print('square, ab = 0 = cd       ', hh([1, 2, 3, 4], sq, {('a', 'b'), ('c', 'd')}), ' (control: different algebra)')
