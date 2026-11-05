"""Round 118 (theorist), T5 / H-023, loop models M4 and M5 of E-284 (workshop/rounds/106/theorist_loopfilt.py).
Proof scheme (see theorist.md): triangular change of generators on kM  =>  kM = B * T(E), E contractible, B = T(X + Theta') a dg subalgebra;
then B -> B/(Theta') is a q.i. by a g-filtration and Lemma A (E-240).  This script checks, on the engine's printed differentials,
(1) every new generator's d equals the claimed expression in the new generators (words multiplied out in the old letters),
(2) the change of generators is invertible (each new generator = old generator + words in other letters of the same (src, tgt, deg, weight)),
(3) g never raises under d (control: a planted g that does raise is caught).
Usage: .venv/bin/python -W ignore workshop/rounds/118/theorist_m45.py"""
import importlib.util, itertools
sp = importlib.util.spec_from_file_location('lf', 'workshop/rounds/106/theorist_loopfilt.py')
lf = importlib.util.module_from_spec(sp); sp.loader.exec_module(lf)
def add(*ds):
    out = {}
    for c, d in ds:
        for w, k in d.items(): out[w] = out.get(w, 0) + c * k
    return {w: k for w, k in out.items() if k}
def mul(x, y): return {a + b: p * q for a, p in x.items() for b, q in y.items()}
def G(L): return {(L,): 1}

def run(name, repl, newgen, claim, gdef, gplant):
    """repl: new generator -> old letter it replaces (or itself); newgen: new generator -> dict word(old letters) -> coef;
    claim: new generator -> dict word(new generators) -> coef = its differential."""
    q = lf.MODELS[name](); m, letters, wt = lf.setup(q); ok = True
    def subst(expr):
        out = {}
        for w, k in expr.items():
            cur = {(): k}
            for x in w: cur = mul(cur, newgen[x])
            out = add((1, out), (1, cur))
        return out
    for x in newgen:
        lhs = {w: k for w, k in m.partial(newgen[x]).items() if k}; rhs = subst(claim[x]); good = lhs == rhs; ok &= good
        print(f'  d {x:5s} = {claim[x]}  [{"ok" if good else "MISMATCH"}]')
    # (2a) homogeneity, (2b) leading term: the replaced old letter occurs once with coefficient +-1, (2c) dependency graph acyclic
    dep = {}
    for x, e in newgen.items():
        metas = {(m.src(w[0]), m.tgt(w[-1]), sum(m.deg(Lt) for Lt in w), sum(wt(Lt) for Lt in w)) for w in e}
        if len(metas) != 1: ok = False; print('  INHOMOGENEOUS', x, metas)
        lead = e.get((repl[x],), 0)
        if abs(lead) != 1: ok = False; print('  NO UNIT LEADING TERM', x)
        dep[repl[x]] = {Lt for w in e for Lt in w if Lt != repl[x]}
    def acyc(n, seen=()):
        if n in seen: return False
        return all(acyc(k, seen + (n,)) for k in dep.get(n, ()) if k in dep)
    cyc = [n for n in dep if not acyc(n)]
    print('  substitution dependency graph acyclic:', not cyc); ok &= not cyc
    for gname, g in (('g', gdef), ('planted g', gplant)):
        raised = [(x, w) for x in claim for w in claim[x] if sum(g.get(y, 0) for y in w) > g.get(x, 0)]
        print(f'  {gname}: d raises g on', raised if raised else 'nothing')
        ok &= (not raised) if gname == 'g' else bool(raised)
    print(' ', name, 'ALL OK' if ok else 'FAILED'); return ok
def Ls(*a): return ('m',) + a
def lin(*terms):   # terms: (coef, tuple of letters)
    return add(*[(c, {w: 1}) for c, w in terms])
def nw(*terms): return lin(*terms)
# ---------- M4: letters ----------
b = Ls(None, 'b', None); c1 = Ls('a1', 'b', None); c2 = Ls('a2', 'b', None); fb = Ls(None, 'f', None)
u1 = Ls('a1', 'f', None); u2 = Ls('a2', 'f', None); t1 = Ls(None, 'f', 'a1'); t2 = Ls(None, 'f', 'a2')
s11 = Ls('a1', 'f', 'a1'); s12 = Ls('a1', 'f', 'a2'); s21 = Ls('a2', 'f', 'a1'); s22 = Ls('a2', 'f', 'a2'); v1 = ('s', 'a1'); v2 = ('s', 'a2')
d_t1 = lin((1, (s11, v1)), (1, (s21, v2)), (-1, (v1, fb)), (1, (b,)))
d_s11 = lin((1, (v1, u1)), (-1, (c1,))); d_s21 = lin((1, (v1, u2)), (-1, (c2,)))
newgen4 = {'p1': d_t1, 'q11': d_s11, 'q21': d_s21, 't1': nw((1, (t1,))), 's11': nw((1, (s11,))), 's21': nw((1, (s21,))),
  'sg': lin((1, (t1,)), (1, (t2,))), 'r22': lin((1, (s22,)), (1, (s21,))), 'sp': lin((1, (s12,)), (1, (s11,)), (-1, (s22,)), (-1, (s21,))),
  'w': lin((1, (v1,)), (1, (v2,))), 'v1': nw((1, (v1,))), 'u1p': lin((1, (u1,)), (-1, (u2,))), 'u2': nw((1, (u2,))), 'fb': nw((1, (fb,)))}
repl4 = {'p1': b, 'q11': c1, 'q21': c2, 't1': t1, 's11': s11, 's21': s21, 'sg': t2, 'r22': s22, 'sp': s12, 'w': v2, 'v1': v1, 'u1p': u1, 'u2': u2, 'fb': fb}
def C(**kw): return kw
claim4 = {'p1': {}, 'q11': {}, 'q21': {}, 't1': {('p1',): 1}, 's11': {('q11',): 1}, 's21': {('q21',): 1},
  'sg': {('sp', 'v1'): 1, ('r22', 'w'): 1, ('w', 'fb'): -1}, 'r22': {('w', 'u2'): 1}, 'sp': {('w', 'u1p'): 1}, 'w': {}, 'v1': {}, 'u1p': {}, 'u2': {},
  'fb': {('u1p', 'v1'): 1, ('u2', 'w'): 1}}
g4 = {'sg': 1, 'fb': 1}; g4bad = {'sg': 1, 'fb': 1, 'w': 1}
print('== M4'); ok4 = run('M4_two_alpha', repl4, newgen4, claim4, g4, g4bad)
# ---------- M5 ----------
e = Ls(None, 'e', None); ca = Ls('a', 'b', None); cm = Ls('m', 'b', None); ua = Ls('a', 'f', None); um = Ls('m', 'f', None)
ta = Ls(None, 'f', 'a'); tm = Ls(None, 'f', 'm'); saa = Ls('a', 'f', 'a'); sam = Ls('a', 'f', 'm'); sma = Ls('m', 'f', 'a'); smm = Ls('m', 'f', 'm')
va = ('s', 'a'); vm = ('s', 'm')
d_ta = lin((1, (saa, va)), (1, (sma, vm)), (-1, (va, fb)), (1, (b,)))
d_saa = lin((1, (va, ua)), (-1, (ca,))); d_sma = lin((1, (va, um)), (-1, (cm,)))
newgen5 = {'pa': d_ta, 'qaa': d_saa, 'qma': d_sma, 'ta': nw((1, (ta,))), 'saa': nw((1, (saa,))), 'sma': nw((1, (sma,))),
  'tmp': lin((1, (tm,)), (-1, (e, ta))), 'smmp': lin((1, (smm,)), (-1, (e, sma))),
  'w': lin((1, (vm,)), (-1, (e, va))),
  'samp': lin((1, (sam,)), (-1, (e, saa)), (1, (smm, e)), (-1, (e, sma, e))),
  'uap': lin((1, (ua,)), (1, (um, e))), 'um': nw((1, (um,))), 'va': nw((1, (va,))), 'eps': nw((1, (e,))), 'fb': nw((1, (fb,)))}
repl5 = {'pa': b, 'qaa': ca, 'qma': cm, 'ta': ta, 'saa': saa, 'sma': sma, 'tmp': tm, 'smmp': smm, 'w': vm, 'samp': sam, 'uap': ua, 'um': um, 'va': va, 'eps': e, 'fb': fb}
claim5 = {'pa': {}, 'qaa': {}, 'qma': {}, 'ta': {('pa',): 1}, 'saa': {('qaa',): 1}, 'sma': {('qma',): 1},
  'tmp': {('samp', 'va'): 1, ('smmp', 'w'): 1, ('w', 'fb'): -1}, 'smmp': {('w', 'um'): 1}, 'samp': {('w', 'uap'): 1}, 'w': {}, 'uap': {}, 'um': {}, 'va': {}, 'eps': {},
  'fb': {('uap', 'va'): 1, ('um', 'w'): 1}}
g5 = {'tmp': 1, 'fb': 1}; g5bad = {'tmp': 1, 'fb': 1, 'w': 1}
print('== M5'); ok5 = run('M5_loop_two_paths', repl5, newgen5, claim5, g5, g5bad)
print('M4', ok4, 'M5', ok5)
