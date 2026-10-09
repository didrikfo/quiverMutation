"""Resolve the cospectral-unresolved LNAs of length n by the F-047 profile (Smith form of g(Phi) per irreducible factor g of the Coxeter
polynomial, plus Smith form of C + C^T): a derived invariant independent of the Coxeter key and of the move orbits.
Soundness check: the profile is constant on every resolved class of the cospectral key; separation check: seed profiles differ.
usage: python workshop/rounds/029/toolsmith_snfresolve.py N [--plan]"""
import sys, collections, time
import sympy
from sympy.matrices.normalforms import smith_normal_form
sys.path.insert(0, "workshop/rounds/027")
import maverick_classes as mc
from quivermutation import coxeterTables as ct, quipuForms as qf, nakayama as nk
T = sympy.symbols('T')
_cache = {}
def profile(n, rl):
    C = sympy.Matrix(ct.lnaCartanMatrix(n, rl))
    k = tuple(C)
    if k in _cache: return _cache[k]
    Phi = -(C.inv().T) * C
    cp = sympy.Poly(Phi.charpoly(T).as_expr(), T)
    snf = lambda M: tuple(sorted(abs(int(s[i, i])) for s in [smith_normal_form(M, domain=sympy.ZZ)] for i in range(n)))
    out = [snf(C + C.T)]
    for g, e in sympy.factor_list(cp.as_expr())[1]:
        G = sympy.zeros(n, n)
        for (d,), c in sympy.Poly(g, T).terms(): G += c * Phi**d
        out.append((str(g), e, snf(G)))
    _cache[k] = tuple(out); return _cache[k]
if __name__ == '__main__':
    n = int(sys.argv[1]); plan = '--plan' in sys.argv
    lab, bad = mc.classes(n)
    print(n, 'unresolved', len(bad), 'distinct unresolved Cartans', len({tuple(map(tuple, ct.lnaCartanMatrix(n, m))) for m in bad}))
    keys = collections.Counter(lab[m][0] for m in lab if lab[m][1] == '?'); print('keys of unresolved', len(keys), dict(keys) if len(keys) < 5 else '')
    if plan:
        t = time.time(); profile(n, bad[0]); print('one profile %.2f s; estimate %.0f s' % (time.time() - t, (time.time() - t) * len({tuple(map(tuple, ct.lnaCartanMatrix(n, m))) for m in bad}))); sys.exit()
    # per cospectral key: resolved classes -> set of profiles, then unresolved matched against them
    reskey = collections.defaultdict(lambda: collections.defaultdict(set)); t0 = time.time()
    keyset = {lab[m][0] for m in bad}
    for m, l in lab.items():
        if l[0] in keyset and l[1] != '?': reskey[l[0]][l[1]].add(profile(n, m))
    for k, d in reskey.items():
        for c, ps in d.items(): print('key', k, 'class', c, 'profiles', len(ps), 'sound' if len(ps) == 1 else 'NOT CONSTANT')
    prof2class = {}
    for k, d in reskey.items():
        for c, ps in d.items():
            for p in ps: prof2class.setdefault((k, p), set()).add(c)
    res = collections.Counter(); out = {}
    for m in bad:
        k = lab[m][0]; cs = prof2class.get((k, profile(n, m)))
        out[m] = cs; res[tuple(sorted(cs)) if cs else None] += 1
    print('result', dict(res)); print('time %.0f s' % (time.time() - t0))
    for m in bad[:20]: print(''.join(map(str, m)), out[m])
