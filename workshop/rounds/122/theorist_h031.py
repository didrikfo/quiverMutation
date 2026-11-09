"""Round 122 (theorist), H-031 as three lemmas: checks of the hypotheses on loop models (formal dg models only).
Engine: workshop/rounds/106/theorist_loopfilt.py (Oppermann 1504.02617v1 Def. 3.2, signs last = -1, E-202/E-284).
Parts (each prints a verdict line; every check has a control that must say NO):
  [LF]   local finiteness: every oriented cycle of the letter quiver of kM has positive weight (=> each (src,tgt,weight)
         cell is finite-dimensional with bounded word length).  Control: a toy letter quiver with a weight-0 cycle.
  [H1]   linear part of d: only Theta letters have a linear term, no linear term is a Theta letter, the pivot letters
         are distinct with unit coefficient (no linear change needed).  Control: M4 (needs a linear change, E-309).
  [PI]   Lemma 2 data: the projection pi: kM -> T(H) killing E = {pivot thetas, their d}, pivot letters solved from
         pi(d theta) = 0; checks pi is a dg map (pi d x = d_H pi x on EVERY letter) and d_H^2 = 0; prints d_H on
         the reduced Theta letters.  Control: a planted sign in one pivot solution must break the dg-map check.
  [PSI]  Lemma 1(c) instance: the parent automorphism a -> a + m e carries M5' (d f = ab, spectators m, e) to M5
         (d f = ab + meb); the induced letter map Psi: kM(M5') -> kM(M5) commutes with d on all letters and sends
         Theta' into the ideal (Theta).  Control: drop one term of Psi(s_am) -> must fail.
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/122/theorist_h031.py"""
import importlib.util, sys
sp = importlib.util.spec_from_file_location('lf', 'workshop/rounds/106/theorist_loopfilt.py')
lf = importlib.util.module_from_spec(sp); sp.loader.exec_module(lf)
from fractions import Fraction as Fr
sp2 = importlib.util.spec_from_file_location('lf2', 'workshop/rounds/106/theorist_loopfilt2.py')
lf2 = importlib.util.module_from_spec(sp2); sp2.loader.exec_module(lf2)

MODELS = dict(lf.MODELS); MODELS.update(lf2.M)
# M5': M1 plus spectator arrows m: i->k, e: k->j (d f = ab only); same weights as M5
MODELS['M5p_M1_plus_spectators'] = lambda: lf.mk({'a': (0, 1, 0), 'b': (1, 0, 0), 'f': (0, 0, 1), 'm': (0, 2, 0), 'e': (2, 1, 0)}, ['a', 'm'],
    {'f': {('a', 'b'): 1}}, {'a': 2, 'b': 2, 'f': 4, 'm': 1, 'e': 1})
# M9: d f = ab + m e c with a SECOND arrow c: j->i (parallel to b, c != b), e: k->j; r_m = e c has no pivot and is not absorbed by
# any a -> a + (path), because its last arrow c is not the last arrow of a pivot term. Three vertices (the engine's lf.run scans range(3)).
MODELS['M9_unabsorbable'] = lambda: lf.mk({'a': (0, 1, 0), 'b': (1, 0, 0), 'f': (0, 0, 1), 'm': (0, 2, 0), 'e': (2, 1, 0), 'c': (1, 0, 0)}, ['a', 'm'],
    {'f': {('a', 'b'): 1, ('m', 'e', 'c'): 1}}, {'a': 2, 'b': 2, 'f': 4, 'm': 1, 'e': 1, 'c': 2})

def nm(L):
    if isinstance(L, str): return L
    if L[0] == 's': return L[1] + '*'
    _, a, p, b = L
    return (a or '') + ('.' if a else '') + p + ('/' + b if b else '')
def show(E):
    return ' + '.join(f'{"" if k == 1 else ("-" if k == -1 else str(k) + " ")}{"·".join(nm(x) for x in w) or "1"}' for w, k in sorted(E.items(), key=lambda t: repr(t[0]))) or '0'
def addto(D, w, c):
    D[w] = D.get(w, 0) + c
    if D[w] == 0: del D[w]
def mul(E, F):
    out = {}
    for u, c in E.items():
        for v, e in F.items(): addto(out, u + v, c * e)
    return out

def build(name):
    q = MODELS[name](); m, letters, wt = lf.setup(q)
    dval = {L: {w: k for w, k in m.partial({(L,): 1}).items() if k} for L in letters}
    theta = {L for L in letters if L[0] == 'm' and L[3] is not None}
    return q, m, letters, wt, dval, theta

def leib(dval, deg, E):   # d on linear combinations of words, applied-first sign rule of the engine
    out = {}
    for w, c in E.items():
        for k, L in enumerate(w):
            s = (-1) ** sum(deg(x) for x in w[k + 1:])
            for u, e in dval.get(L, {}).items(): addto(out, w[:k] + u + w[k + 1:], c * s * e)
    return out

# ---------------- [LF]
def lf_check(verts, edges):
    INF = float('inf'); D = {(u, v): INF for u in verts for v in verts}
    for s, t, w in edges: D[(s, t)] = min(D[(s, t)], w)
    for k in verts:
        for u in verts:
            for v in verts:
                if D[(u, k)] + D[(k, v)] < D[(u, v)]: D[(u, v)] = D[(u, k)] + D[(k, v)]
    mins = {v: D[(v, v)] for v in verts}
    return all(x > 0 for x in mins.values()), mins

# ---------------- [H1] + [PI]
def lin_part(dval, L): return {w[0]: k for w, k in dval[L].items() if len(w) == 1}

def h1(letters, dval, theta):
    src = [L for L in letters if lin_part(dval, L)]
    nontheta_src = [nm(L) for L in src if L not in theta]
    theta_tgt = [nm(x) for L in src for x in lin_part(dval, L) if x in theta]
    piv = {}; bad = []
    for L in sorted(src, key=repr):
        lp = lin_part(dval, L)
        cand = [x for x, k in lp.items() if abs(k) == 1 and x not in theta]
        x = next((x for x in cand if x not in piv.values() and all(x not in lin_part(dval, M) for M in src if M != L)), None)
        if x is None: bad.append(nm(L))
        else: piv[L] = x
    return dict(nontheta_src=nontheta_src, theta_tgt=theta_tgt, unpivoted=bad), piv

def projection(q, m, letters, dval, theta, piv, plant=False):
    deg = m.deg; Epiv = set(piv); pivlet = {x: L for L, x in piv.items()}
    sol = {}   # pivot letter x -> expression with pi(d theta) = 0
    for L, x in piv.items():
        c = dval[L][(x,)]; rest = {w: k for w, k in dval[L].items() if w != (x,)}
        sol[x] = {w: Fr(-k, c) for w, k in rest.items()}
        if plant and not plant_done[0]: sol[x] = {w: -k for w, k in sol[x].items()}; plant_done[0] = True
    def pi(E, depth=0):
        assert depth < 60, 'pivot substitution does not terminate'
        out = {}; again = False
        for w, k in E.items():
            if any(x in Epiv for x in w): continue
            cur = {(): Fr(k)}
            for x in w:
                if x in sol: cur = mul(cur, sol[x]); again = True
                else: cur = mul(cur, {(x,): 1})
            for u, e in cur.items(): addto(out, u, e)
        out = {w: k for w, k in out.items() if not any(x in Epiv for x in w)}
        return pi(out, depth + 1) if again else out
    H = [L for L in letters if L not in Epiv and L not in sol]
    dH = {L: pi(dval[L]) for L in H}
    bad_dg = [nm(x) for x in letters if (pi(leib(dval, deg, {(x,): 1})) if True else None) != leib(dH, deg, pi({(x,): 1}))]
    bad_d2 = [nm(L) for L in H if leib(dH, deg, dH[L])]
    return H, dH, bad_dg, bad_d2
plant_done = [False]

def part_models(names):
    for name in names:
        q, m, letters, wt, dval, theta = build(name)
        verts = sorted({m.src(L) for L in letters} | {m.tgt(L) for L in letters})
        ok, mins = lf_check(verts, [(m.src(L), m.tgt(L), wt(L)) for L in letters])
        r, piv = h1(letters, dval, theta)
        h1ok = not r['nontheta_src'] and not r['theta_tgt'] and not r['unpivoted']
        print(f'== {name}: {len(letters)} letters, {len(theta)} Theta')
        print(f'  [LF] min closed-walk weight per vertex {mins} -> {"YES" if ok else "NO"}')
        print(f'  [H1] {r} pivots {[(nm(L), nm(x)) for L, x in piv.items()]} -> {"YES" if h1ok else "NO"}')
        if not h1ok: continue
        H, dH, bad_dg, bad_d2 = projection(q, m, letters, dval, theta, piv)
        print(f'  [PI] H = {len(H)} letters; pi dg map fails on {bad_dg}; d_H^2 != 0 on {bad_d2} -> {"YES" if not bad_dg and not bad_d2 else "NO"}')
        for L in H:
            if L in theta or dH[L]: print(f'     d_H {nm(L):8s}{"(Theta)" if L in theta else "       "} = {show(dH[L])}')
        if name == 'M5_loop_two_paths':
            plant_done[0] = False
            _, _, bdg, _ = projection(q, m, letters, dval, theta, piv, plant=True)
            print(f'  [PI control] planted sign in one pivot solution: dg-map failures {len(bdg)} -> {"NO (caught)" if bdg else "YES (NOT caught)"}')

# ---------------- [PSI]
def part_psi():
    _, mp, Lp, _, dp, thp = build('M5p_M1_plus_spectators')
    q5, m5, L5, _, d5, th5 = build('M5_loop_two_paths')
    def X(a, p, b): return ('m', a, p, b)
    e = X(None, 'e', None); va, vm = ('s', 'a'), ('s', 'm')
    def P(*terms): return {tuple(w): k for k, w in terms}
    def psi_map(drop=False):
        S = {X(None, 'b', None): P((1, [X(None, 'b', None)])), X('a', 'b', None): P((1, [X('a', 'b', None)]), (1, [X('m', 'b', None), e])),
             X('m', 'b', None): P((1, [X('m', 'b', None)])), X(None, 'f', None): P((1, [X(None, 'f', None)])),
             X('a', 'f', None): P((1, [X('a', 'f', None)]), (1, [X('m', 'f', None), e])), X('m', 'f', None): P((1, [X('m', 'f', None)])),
             X(None, 'f', 'a'): P((1, [X(None, 'f', 'a')])), X(None, 'f', 'm'): P((1, [X(None, 'f', 'm')]), (-1, [e, X(None, 'f', 'a')])),
             X('a', 'f', 'a'): P((1, [X('a', 'f', 'a')]), (1, [X('m', 'f', 'a'), e])), X('m', 'f', 'a'): P((1, [X('m', 'f', 'a')])),
             X('m', 'f', 'm'): P((1, [X('m', 'f', 'm')]), (-1, [e, X('m', 'f', 'a')])),
             X('a', 'f', 'm'): P((1, [X('a', 'f', 'm')]), (-1, [e, X('a', 'f', 'a')]), (1, [X('m', 'f', 'm'), e]), (-1, [e, X('m', 'f', 'a'), e])),
             e: P((1, [e])), va: P((1, [va])), vm: P((1, [vm]), (-1, [e, va]))}
        if drop: S[X('a', 'f', 'm')].pop((X('m', 'f', 'm'), e))
        return S
    for drop in (False, True):
        S = psi_map(drop)
        missing = [nm(L) for L in Lp if L not in S]
        def Psi(E):
            out = {}
            for w, k in E.items():
                cur = {(): k}
                for x in w: cur = mul(cur, S[x])
                for u, c in cur.items(): addto(out, u, c)
            return out
        bad = [nm(L) for L in Lp if L in S and Psi(dp[L]) != leib(d5, m5.deg, S[L])]
        th_ok = all(all(any(x in th5 for x in w) for w in S[L]) for L in thp)
        tag = 'control (dropped smm.e in Psi(s_am))' if drop else 'Psi'
        print(f'[PSI] {tag}: letters unmapped {missing}; d Psi != Psi d on {bad}; Psi(Theta\') in (Theta): {th_ok} -> {"YES" if not bad and not missing else "NO"}')
        if not drop:   # inverse: same formulas with e -> -e applied to M5 letters gives M5' elements; check Psi o Psi^-1 = id
            Sinv = {L: {w: k * (-1) ** sum(1 for x in w if x == e) for w, k in S[L].items()} for L in S}
            # Psi^-1(L) = L - (corrections): solve triangularly by iteration
            inv = {}
            for L in S:
                cur = {(L,): 1}
                for _ in range(6):
                    img = Psi(cur); diff = {w: -k for w, k in img.items()}; addto(diff, (L,), 1)
                    if not diff: break
                    for w, k in diff.items(): addto(cur, w, k)
                inv[L] = cur
            ok_inv = all(Psi(inv[L]) == {(L,): 1} for L in S)
            th_inv = all(all(any(x in thp for x in w) for w in inv[L]) for L in th5)
            print(f'[PSI] Psi invertible (triangular inverse found, Psi(Psi^-1 L) = L on all {len(S)} letters): {ok_inv}; Psi^-1(Theta) in (Theta\'): {th_inv}')

if __name__ == '__main__':
    ok, mins = lf_check([0, 1], [(0, 1, 2), (1, 0, -2)])
    print(f'[LF control] toy quiver with a weight-0 cycle -> {"YES (NOT caught)" if ok else "NO (caught)"}')
    part_models(['M1_loop1_dphi=ba', 'M4_two_alpha', 'M5_loop_two_paths', 'M5p_M1_plus_spectators', 'M9_unabsorbable', 'M3b_loop2_dfnonzero'])
    part_psi()

# ---------------- [SH] Lemma 3 shape on (T(H), d_H): search g: X -> {0..G}; theta's g = g of its leading word
def shape(m, H, dH, theta, g):
    X = [L for L in H if L not in theta]; Th = [L for L in H if L in theta]
    G = dict(g)
    for x in X:
        if any(sum(G[y] for y in w) >= G[x] for w in dH[x] if not set(w) & theta): return None, f'd0 != 0 on {nm(x)}'
    lead = {}; todo = list(Th)
    for _ in range(len(Th) + 1):
        for t in list(todo):
            if any(y in theta and y not in G for w in dH[t] for y in w): continue
            gs = sorted(((sum(G[y] for y in w), w) for w in dH[t]), key=lambda z: -z[0])
            top = [w for s, w in gs if s == gs[0][0]]
            if len(top) != 1 or len(top[0]) != 2 or any(y in theta for y in top[0]): return None, f'leading part of d {nm(t)} not one quadratic X-monomial'
            G[t] = gs[0][0]; lead[t] = top[0]; todo.remove(t)
    if todo: return None, 'theta g cyclic'
    for x in X:
        if any(sum(G[y] for y in w) >= G[x] for w in dH[x]): return None, f'd0 != 0 on {nm(x)} (Theta term)'
    X1 = {w[0] for w in lead.values()}; Y = {w[1] for w in lead.values()}
    if X1 & Y: return None, 'X1, Y not disjoint'
    pairs = {(x, y) for x in X1 for y in Y if m.tgt(x) == m.src(y)}
    if pairs != set(lead.values()) or len(set(lead.values())) != len(lead): return None, f'composable pairs {len(pairs)} != thetas {len(lead)}'
    return (G, lead), 'ok'

def search(m, H, dH, theta, Gmax=3):
    import itertools
    X = [L for L in H if L not in theta]; nX = [x for x in X if dH[x]]   # letters with d = 0 can keep g = 0 or 1
    found = None; tried = 0
    for vals in itertools.product(range(Gmax + 1), repeat=len(X)):
        tried += 1
        r, why = shape(m, H, dH, theta, dict(zip(X, vals)))
        if r: found = r; break
    return found, tried

SHAPE = ['M1_loop1_dphi=ba', 'M2_loop1_dphi=0', 'M5_loop_two_paths', 'M5p_M1_plus_spectators', 'M9_unabsorbable', 'M3b_loop2_dfnonzero',
         'M3_loop2_junction', 'M7_two_loops', 'M7b_two_loops_both', 'M8_loop3_A', 'M8_loop3_B']
def part_shape():
    out = {}
    for name in SHAPE:
        q, m, letters, wt, dval, theta = build(name); r1, piv = h1(letters, dval, theta)
        if r1['unpivoted'] or r1['nontheta_src'] or r1['theta_tgt']: print(f'[SH] {name}: H1 fails {r1} -> not reached'); continue
        H, dH, bdg, bd2 = projection(q, m, letters, dval, theta, piv)
        if bdg or bd2: print(f'[SH] {name}: projection not a dg map -> not reached'); continue
        found, tried = search(m, H, dH, theta, 2 if len(H) > 10 else 3)
        if found:
            G, lead = found
            print(f'[SH] {name}: Lemma A shape found after {tried} g-vectors: g = { {nm(k): v for k, v in G.items() if v} }; leading words {[(nm(t), nm(w[0]) + "·" + nm(w[1])) for t, w in lead.items()]}')
            # Bergman ambiguity check: leading words x·y, an overlap needs y in X1 (excluded by X1 cap Y = 0) -> none
            if not lead: print('     (no Theta left in H: Lemma 3 vacuous)'); out[name] = found; continue
            bad = {k: v for k, v in G.items() if k not in theta}; x0 = max(bad, key=lambda k: bad[k]); bad[x0] = 0
            r, why = shape(m, H, dH, theta, bad)
            print(f'     control: planted g({nm(x0)}) = 0 -> {"YES (NOT caught)" if r else "NO (caught): " + why}')
        else:
            print(f'[SH] {name}: no g in the box gives Lemma A shape on the raw letters of H ({tried} tried) -> NO' + (' (no Theta left: Lemma 3 vacuous)' if not any(L in theta for L in H) else ''))
        out[name] = found
    return out

def part_engine(W):
    for name in ['M5p_M1_plus_spectators', 'M9_unabsorbable']:
        lf.run(name, MODELS[name](), W)

def part_nearmiss():
    # E-284 near-miss: the closed ideal I + (a*) of M1 has H != 0 in 5 cells (weight <= 9); treat a* as an extra Theta letter
    q, m, letters, wt, dval, theta = build('M1_loop1_dphi=ba'); th2 = set(theta) | {('s', 'a')}
    _, piv = h1(letters, dval, th2); H, dH, bd, b2 = projection(q, m, letters, dval, th2, piv)
    found, tried = search(m, H, dH, th2, 3)
    print(f'[SH control] M1 with Theta + a* (E-284 near-miss, H(I) != 0 there): pi dg {not bd}; shape found: {bool(found)} ({tried} g tried) -> {"YES (NOT caught)" if found else "NO (caught)"}')

if __name__ == '__main__':
    part_shape(); part_nearmiss()
    if '--engine' in sys.argv: part_engine(int(sys.argv[sys.argv.index('--engine') + 1]))
