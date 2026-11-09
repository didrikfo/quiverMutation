"""Round 122 (scholar), H-031 / T5: test, on the raw kM of the loop models (workshop/rounds/106/theorist_loopfilt.py,
Oppermann 1504.02617v1 Def. 3.2, signs 'last' = -1), the hypotheses of the NON-COMPLETED semifree destabilisation argument
of Etnyre-Ng-Sabloff, arXiv:math/0101145v3 (Def. 3.8 elementary/tame, Def. 3.9 stabilisation, Cor. 3.11 homology unchanged by
stabilisation, Sec. 6.3 'Move III' bootstrap g_i, Lemmas 6.2-6.4; numbers by workshop/rounds/110/scholar_thmnumbers.py) and Chekanov (math/9709233).  Their hypotheses, transcribed to T_K(V), K = k^n:
  (T) triangularity: the generators can be ordered so that d(g) lies in the subalgebra of earlier generators
      (ENS: heights, Lemma 6.3 'd b_i contains only b_j with j > i' and 'd a_i in A_{i-1}');  checked as: the digraph
      g -> letters occurring in d(g) is acyclic;
  (P) pair: d(theta) = c b + v with c != 0, b a generator, v not involving theta or b and no letter of v having b
      below it (so the order low < b < theta < high exists; ENS eq. height-order).
Under (T)+(P), ENS's proof gives a finite tame isomorphism A ~ S(A') (no completion); A' is identified with the dg quotient
A/(theta, d theta), i.e. the substitution theta -> 0, b -> -v/c; we compute it and iterate while some pair satisfies (P).
Checks: d^2 = 0 on every intermediate algebra (Leibniz, applied-first words, sign (-1)^{deg of later letters} as the engine);
remaining letters and Theta letters at the end.  Power: a planted cycle (T fails) and a planted pair whose v contains a letter
lying over b (P fails) must be refused.
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/122/scholar_ens.py"""
import importlib.util
from fractions import Fraction as F
sp = importlib.util.spec_from_file_location('lf', 'workshop/rounds/106/theorist_loopfilt.py')
lf = importlib.util.module_from_spec(sp); sp.loader.exec_module(lf)

def name(L):
    if isinstance(L, str): return L
    if L[0] == 's': return L[1] + '*'
    _, a, p, b = L
    return (a or '') + p + ('/' + b if b else '')

def add(D, w, c):
    D[w] = D.get(w, 0) + c
    if D[w] == 0: del D[w]

def dword(w, d, deg):
    out = {}
    for k, x in enumerate(w):
        s = (-1) ** sum(deg[y] for y in w[k + 1:])
        for u, c in d[x].items(): add(out, w[:k] + u + w[k + 1:], s * c)
    return out

def d2_ok(gens, d, deg):
    for g in gens:
        tot = {}
        for w, c in d[g].items():
            for u, e in dword(w, d, deg).items(): add(tot, u, c * e)
        if tot: return False
    return True

def below(gens, d):
    """letters strictly below g (reachable through d); None if the digraph has a cycle."""
    succ = {g: {x for w in d[g] for x in w} for g in gens}
    memo, state = {}, {}
    def rec(g):
        if state.get(g) == 1: raise ValueError
        if g in memo: return memo[g]
        state[g] = 1; r = set()
        for x in succ[g]: r |= {x} | rec(x)
        state[g] = 2; memo[g] = r; return r
    try: return {g: rec(g) for g in gens}
    except ValueError: return None

def find_pair(gens, d, B, prefer):
    for th in sorted(gens, key=lambda g: (g not in prefer, name(g))):
        for w, c in d[th].items():
            if len(w) != 1 or w[0] == th: continue
            b = w[0]; v = {u: e for u, e in d[th].items() if u != w}
            lv = {x for u in v for x in u}
            if th in lv or b in lv: continue
            if any(b in B[x] for x in lv): continue
            return th, b, c, v
    return None

def subst(E, rep):
    out = {}
    for w, c in E.items():
        acc = {(): F(c)}
        for x in w:
            img = rep.get(x, {(x,): 1}); new = {}
            for u, e in acc.items():
                for z, f in img.items(): add(new, u + z, e * f)
            acc = new
        for u, e in acc.items(): add(out, u, e)
    return out

def destabilise(gens, d, deg, theta, label):
    gens = list(gens); steps = []
    while True:
        B = below(gens, d)
        if B is None: return gens, d, steps, 'NOT TRIANGULAR'
        p = find_pair(gens, d, B, theta)
        if p is None: break
        th, b, c, v = p
        rep = {th: {}, b: {u: F(-e) / c for u, e in v.items()}}
        gens = [g for g in gens if g not in (th, b)]
        d = {g: subst(d[g], rep) for g in gens}
        ok = d2_ok(gens, d, deg)
        steps.append(f'{name(th)} -> {name(b)} (v: {len(v)} words) d^2=0:{ok}')
        if not ok: return gens, d, steps, 'd^2 != 0'
    return gens, d, steps, 'ok'

def model(mod):
    q = lf.MODELS[mod](); m, letters, wt = lf.setup(q)
    d = {L: {w: k for w, k in m.partial({(L,): 1}).items() if k} for L in letters}
    deg = {L: m.deg(L) for L in letters}
    theta = {L for L in letters if L[0] == 'm' and L[3] is not None}
    return letters, d, deg, theta

if __name__ == '__main__':
    for mod in ['M1_loop1_dphi=ba', 'M2_loop1_dphi=0', 'M3_loop2_junction', 'M3b_loop2_dfnonzero', 'M4_two_alpha',
                'M5_loop_two_paths', 'M6_noloop_control']:
        letters, d, deg, theta = model(mod)
        tri = below(letters, d) is not None
        g2, d2, steps, st = destabilise(letters, d, deg, theta, mod)
        rest_th = sorted(name(x) for x in g2 if x in theta)
        quad = all(len(w) >= 2 for g in g2 for w in d2[g])
        print(f'== {mod}: {len(letters)} letters; d^2=0 on kM: {d2_ok(letters, d, deg)}; (T) triangular: {tri}; status {st}')
        for s in steps: print('    cancel', s)
        print(f'    left {len(g2)} letters {sorted(name(x) for x in g2)}; Theta left: {rest_th}; d of the rest has no linear term: {quad}')
        if rest_th: print('    d on Theta left:', {name(x): [([name(y) for y in w], str(c)) for w, c in d2[x].items()] for x in g2 if x in theta})
    # power 1: planted cycle -- in M1 add the word (f/a) to d(a*) formally (breaks d^2 but tests (T) only)
    letters, d, deg, theta = model('M1_loop1_dphi=ba')
    fa = [x for x in letters if name(x) == 'f/a'][0]; ast = [x for x in letters if name(x) == 'a*'][0]
    dp = dict(d); dp[ast] = {(fa,): 1}
    print('== PLANT cycle (d a* := f/a in M1): (T) triangular:', below(letters, dp) is not None, '-> refused' if below(letters, dp) is None else '-> MISSED')
    # power 2: planted (P) failure -- d(theta) = b + x with b below x
    gens = ['b', 'x', 't', 'y', 'z']; dd = {'b': {}, 'y': {}, 'z': {}, 'x': {('b', 'y'): 1}, 't': {('b',): 1, ('x', 'z'): 1}}
    fp = find_pair(gens, dd, below(gens, dd), set())
    print('== PLANT pair (d t = b + x z, b below x): pair found:', fp, '-> refused' if fp is None else '-> MISSED')
    # control YES where YES is true: a stabilisation S(A) of ENS Def. 3.6, d e1 = e2, plus the E-309 substituted M5 (claim5)
    import sys, io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sys.path.insert(0, 'workshop/rounds/118'); import theorist_m45 as T
    gens = list(T.claim5); deg5 = {g: 0 for g in gens}
    B = below(gens, T.claim5)
    g2, _, steps, st = destabilise(gens, T.claim5, deg5, {'ta', 'saa', 'sma', 'tmp', 'smmp', 'samp'}, 'claim5')
    print(f'== CONTROL E-309 substituted M5: triangular {B is not None}; cancelled {len(steps)} pairs {[s.split(" (")[0] for s in steps]}; left {len(g2)} {sorted(g2)}')
    g4 = list(T.claim4); g4r, _, st4, _ = destabilise(g4, T.claim4, {g: 0 for g in g4}, set(), 'claim4')
    print(f'== CONTROL E-309 substituted M4: triangular {below(g4, T.claim4) is not None}; cancelled {len(st4)} pairs; left {len(g4r)} {sorted(g4r)}')
