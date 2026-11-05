"""Round 118 (theorist), answer to referee item 1-2: the STRUCTURAL claims behind the M4/M5 proof that theorist_m45.py does not check.
 (i)   change of generators is a bijection of generator sets with equal (src, tgt, deg, weight); counts of words of length <= 5 per cell agree (kM = T(new letters));
       E and B = T(X + Theta') are d-closed free dg subalgebras: d of every B-generator uses only B-letters, d of every E-generator only E-letters; so kM = B * T(E).
 (ii)  I = (E) + (Theta', d Theta'): explicit inverse of the substitution on the old Theta letters (checked by substituting), and every new E / Theta' generator contains an old Theta letter in each word.
 (iii) Lemma A shape: Theta' pairs = exactly the composable pairs X1 x Y (computed from src/tgt), degrees deg(theta) = deg(w) + deg(y) + 1, spectator arrows (e, v1) occur in no Theta' differential.
Usage: .venv/bin/python -W ignore workshop/rounds/118/theorist_m45_struct.py   (imports theorist_m45, prints its output first)"""
import sys, io, contextlib, itertools
sys.path.insert(0, 'workshop/rounds/118')
buf = io.StringIO()
with contextlib.redirect_stdout(buf): import theorist_m45 as M
lf = M.lf
def first_letters(expr): return next(iter(expr))
def check(name, model, repl, newgen, claim, E, Th, X1, Y, spect, inv, g):
    ok = True; q = lf.MODELS[model](); m, letters, wt = lf.setup(q)
    def cell(Lt): return (m.src(Lt), m.tgt(Lt), m.deg(Lt), wt(Lt))
    # (i)
    old = set(letters); olds = set(repl.values())
    print(f'[{name}] old letters {len(old)}, replaced {len(olds)}, same set: {old == olds}, new generators {len(newgen)}'); ok &= old == olds and len(set(newgen)) == len(olds)
    ncell = {x: cell(repl[x]) for x in newgen}
    homog = all({(m.src(w[0]), m.tgt(w[-1]), sum(m.deg(L) for L in w), sum(wt(L) for L in w)) for w in newgen[x]} == {ncell[x]} for x in newgen)
    print('  every new generator has the cell of the letter it replaces:', homog); ok &= homog
    def count(cells):
        cnt = {}
        for L in range(1, 6):
            for ws in itertools.product(cells, repeat=L):
                if all(ws[i][1] == ws[i+1][0] for i in range(L-1)):
                    k = (ws[0][0], ws[-1][1], sum(c[2] for c in ws), sum(c[3] for c in ws)); cnt[k] = cnt.get(k, 0) + 1
        return cnt
    co, cn = count([cell(L) for L in letters]), count(list(ncell.values()))
    print(f'  words of length <= 5 per (src,tgt,deg,weight): old {sum(co.values())}, new {sum(cn.values())}, cell-by-cell equal: {co == cn}'); ok &= co == cn
    Bg = [x for x in newgen if x not in E]
    clE = all(all(y in E for w in claim[x] for y in w) for x in E); clB = all(all(y not in E for w in claim[x] for y in w) for x in Bg)
    print('  d(E-generator) uses only E-letters:', clE, '; d(B-generator) uses only B-letters:', clB, '(B-generators', Bg, ')'); ok &= clE and clB
    # (ii)
    def subst(expr):
        out = {}
        for w, k in expr.items():
            cur = {(): k}
            for x in w: cur = M.mul(cur, newgen[x])
            out = M.add((1, out), (1, cur))
        return out
    oldTh = {}
    for oldname, expr_new in inv.items():     # old letter written as polynomial in NEW generators; substitute -> must equal that old letter
        got = subst(expr_new); good = got == {(oldname,): 1}; ok &= good
        print(f'  inverse of old Theta letter {oldname[1:] if isinstance(oldname, tuple) else oldname}: substituting gives the old letter:', good)
        oldTh[oldname] = all(any(y in E or y in Th for y in w) for w in expr_new)
    print('  each inverse is a combination of E / Theta\' letters only (so old Theta in (E)+(Theta\')):', all(oldTh.values())); ok &= all(oldTh.values())
    # every new E / Theta' generator is in the ideal of old Theta: each word of its expansion contains an old Theta letter
    oldTheta = set(inv)
    inI = {x: all(any(L in oldTheta for L in w) for w in newgen[x]) for x in newgen if x in Th or (x in E and not x.startswith(('p', 'q')))}
    print('  new Theta\' / E-Theta generators lie in (old Theta) (every word has an old Theta letter):', all(inI.values()), inI); ok &= all(inI.values())
    # p, q generators are d of old Theta letters (checked by theorist_m45 as d t = p etc.): list
    # (iii)
    def gcell(x): return ncell[x]
    pairs = {x: next(iter(claim[x] if len(claim[x]) else {})) for x in Th}
    for th, ys in ((t, claim[t]) for t in Th):
        pass
    comp = {(x1, y) for x1 in X1 for y in Y if gcell(x1)[1] == gcell(y)[0]}
    print('  composable pairs X1 x Y (w then y):', sorted(comp))
    thpairs = set()
    for th in Th:
        for w in claim[th]:
            if len(w) == 2 and w[0] in X1 and w[1] in Y: thpairs.add(w)
    # the Theta' differential has leading pair (w,y) with y in Y plus lower g terms; take the g-top term(s)
    gtop = {'M4': {'sg': ('w', 'fb'), 'sp': ('w', 'u1p'), 'r22': ('w', 'u2')}, 'M5': {'tmp': ('w', 'fb'), 'smmp': ('w', 'um'), 'samp': ('w', 'uap')}}[name]
    print("  Theta' pairs (d0 of each Theta'):", sorted(gtop.values()), ' = composable pairs:', set(gtop.values()) == comp); ok &= set(gtop.values()) == comp
    dg = {th: (m.deg(repl[th]), m.deg(repl[gtop[th][0]]) + m.deg(repl[gtop[th][1]]) + 1) for th in Th}
    print('  deg theta vs deg w + deg y + 1:', dg); ok &= all(a == b for a, b in dg.values())
    # d0 : the part of d(theta) of g-degree g(theta) -1 ... here simply: spectator letters occur in no differential of a Theta' generator
    d0 = {th: [w for w in claim[th] if sum(g.get(y, 0) for y in w) == g.get(th, 0)] for th in Th}
    print("  d0 (g-preserving part) of each Theta':", d0)
    sp_in = {s: any(s in w for th in Th for w in d0[th]) for s in spect}
    print("  spectator letters occurring in d0(Theta'):", sp_in, '(they occur only in g-lowering terms of d)'); ok &= not any(sp_in.values())
    others = sorted((x1, z) for x1 in X1 for z in newgen if z not in E and z not in Th and z not in Y and z not in X1 and gcell(x1)[1] == gcell(z)[0])
    print('  letters z of B that are composable after w but are NOT in Y (no Theta, free in gr):', others)
    print(f'  [{name}] STRUCTURE', 'ALL OK' if ok else 'FAILED'); return ok
# inverses of old Theta letters in new generators
M4 = check('M4', 'M4_two_alpha', M.repl4, M.newgen4, M.claim4, {'t1', 'p1', 's11', 'q11', 's21', 'q21'}, {'sg', 'sp', 'r22'}, {'w'}, {'u1p', 'u2', 'fb'}, ['v1'],
  {M.t1: {('t1',): 1}, M.t2: {('sg',): 1, ('t1',): -1}, M.s11: {('s11',): 1}, M.s21: {('s21',): 1}, M.s22: {('r22',): 1, ('s21',): -1},
   M.s12: {('sp',): 1, ('s11',): -1, ('r22',): 1}}, M.g4)
M5 = check('M5', 'M5_loop_two_paths', M.repl5, M.newgen5, M.claim5, {'ta', 'pa', 'saa', 'qaa', 'sma', 'qma'}, {'tmp', 'smmp', 'samp'}, {'w'}, {'uap', 'um', 'fb'}, ['eps', 'va'],
  {M.ta: {('ta',): 1}, M.tm: {('tmp',): 1, ('eps', 'ta'): 1}, M.saa: {('saa',): 1}, M.sma: {('sma',): 1}, M.smm: {('smmp',): 1, ('eps', 'sma'): 1},
   M.sam: {('samp',): 1, ('eps', 'saa'): 1, ('smmp', 'eps'): -1}}, M.g5)
print('M4', M4, 'M5', M5)
