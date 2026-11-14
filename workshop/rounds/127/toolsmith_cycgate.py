"""Round 127 (toolsmith), T5 / H-023: gate-admitted steps on CYCLIC monomial parents, and the (LF), (H1w), (H3), H(I) = 0
checks of E-326 on them.

A loop letter of degree >= 1 at the mutation vertex i forces i onto an oriented cycle (E-279 (1)), so D(i) is cyclic and the
step is never admitted: a gate-admitted *loop-at-i* parent does not exist.  What can exist is a parent whose quiver has an
oriented cycle (with loop letters of degree >= 1 at the cycle's vertices) OUTSIDE D(i), e.g. upstream of i.  This script

  1. enumerates connected quivers on n vertices (no loops, no parallel arrows, at least one oriented cycle), monomial zero
     relations of length 2..LREL (antichain), finite-dimensional, one per isomorphism class;
  2. at every vertex i with an out-arrow and D(i) acyclic, decides the gate EXACTLY (every nonzero path ending at i, not only
     vertex-simple ones) and, separately, the library gate `mutation.mutationIsPossibleAtVertex` and the library J test
     `mutation.isTiltingAtVertex` (both enumerate vertex-simple paths only, arrowPaths.allPathsBetween);
  3. builds the Construction-2.2 minimal model (generic monomial Tate builder, cells (s, t, path-length weight)), checks
     H = A, and runs E-326's analysis (126/theorist_cycparent.analyse: kM at i, d^2, (LF), (C), (H1w), (H3) shape, H(I))
     on every admitted instance, on A (= F step, E-252) and on A^op (= R step).

For a monomial algebra, J = 0 (Ladkani 2.3(c)) is the same as the exact gate: the images p.beta of distinct nonzero paths
are distinct paths, so the map is injective iff no nonzero p is killed by every beta.  Both are therefore one test here.
Controls: [MM] dropped d-term must fail H = A; [GATE] a planted relation that kills c.d on the n = 4 seed must be refused;
[SIMPLE] a seed where the only gate-violating path is non-simple (library admits, exact gate refuses).
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/127/toolsmith_cycgate.py seed
                   .venv/bin/python -W ignore workshop/rounds/127/toolsmith_cycgate.py sweep N MAXARROWS LREL [--plan]"""
import sys, itertools, importlib.util, time
sys.path.insert(0, '.'); sys.path.insert(0, 'tests')
def load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
CP = load('cycp', 'workshop/rounds/126/theorist_cycparent.py'); P = CP.P; h = CP.h; opp = CP.opp
from quivermutation import mutation
from helpers import path_algebra

# ---------------------------------------------------------------- monomial algebra combinatorics
def nz(path, rels):
    """path = tuple of arrow names; nonzero iff no relation occurs as a contiguous subpath"""
    for r in rels:
        L = len(r)
        for o in range(len(path) - L + 1):
            if path[o:o + L] == r: return False
    return True

def nonzero_paths(arrows, rels, maxlen=40):
    """all nonzero paths of length >= 1; None if some nonzero path is longer than maxlen (infinite-dimensional)"""
    out = []; frontier = [(a,) for a in arrows]
    while frontier:
        out += frontier; nxt = []
        for p in frontier:
            t = arrows[p[-1]][1]
            for b, (s, _) in arrows.items():
                if s == t and nz(p + (b,), rels): nxt.append(p + (b,))
        if nxt and len(nxt[0]) > maxlen: return None
        frontier = nxt
    return out

def reach(n, arrows, i):
    seen, st = {i}, [i]
    while st:
        v = st.pop()
        for s, t in arrows.values():
            if s == v and t not in seen: seen.add(t); st.append(t)
    return seen

def has_cycle(verts, arrows):
    E = [(s, t) for s, t in arrows.values() if s in verts and t in verts]
    indeg = {v: 0 for v in verts}
    for s, t in E: indeg[t] += 1
    st = [v for v in verts if indeg[v] == 0]; k = 0
    while st:
        v = st.pop(); k += 1
        for s, t in E:
            if s == v:
                indeg[t] -= 1
                if indeg[t] == 0: st.append(t)
    return k < len(verts)

def exact_gate(arrows, rels, paths, i):
    """the gate read over ALL nonzero paths ending at i (from a vertex != i); also = J = 0 for monomial algebras"""
    outs = [b for b, (s, _) in arrows.items() if s == i]
    if not outs: return False, None
    for p in paths:
        if arrows[p[-1]][1] == i and arrows[p[0]][0] != i and not any(nz(p + (b,), rels) for b in outs):
            return False, p
    return True, None

def library(n, arrows, rels, i):
    names = sorted(arrows)
    vp = lambda path: [arrows[path[0]][0]] + [arrows[a][1] for a in path]
    A = path_algebra([list(arrows[a]) for a in names], [[vp(r)] for r in rels])
    return mutation.mutationIsPossibleAtVertex(A, i), mutation.isTiltingAtVertex(A, i)

def opposite(arrows, rels):
    return {a: (t, s) for a, (s, t) in arrows.items()}, [tuple(reversed(r)) for r in rels]

def gldim_finite(arrows, rels):
    """Anick chains of a monomial algebra: an (m+1)-chain is c_m.s where t_m.s ends with a relation that overlaps t_m
    (starts inside t_m) and t_m.s minus its last arrow contains no relation; the new tail is s.  Extensions depend on the
    tail only, so gl.dim < inf iff the tail graph reachable from the relation tails has no cycle."""
    LM = max(len(r) for r in rels)
    def exts(t):
        out = []
        def rec(sfx):
            if sfx:
                w = t + sfx
                if not nz(w[:-1], rels): return
                if any(w[len(w) - len(r):] == r and len(w) - len(r) < len(t) for r in rels): out.append(sfx); return
            if len(sfx) >= LM - 1: return
            v = arrows[(sfx or t)[-1]][1]
            for b, (s0, _) in sorted(arrows.items()):
                if s0 == v: rec(sfx + (b,))
        rec(()); return out
    state = {}
    def dfs(t):
        state[t] = 1
        for s2 in exts(t):
            if state.get(s2) == 1: return False
            if s2 not in state and not dfs(s2): return False
        state[t] = 2; return True
    return all(dfs(r[1:]) for r in rels if r[1:] not in state)

# ---------------------------------------------------------------- generic monomial Construction-2.2 builder
def minimal_model(n, arrows, rels, WMAX=30, cap=60, plant=False):
    """Tate-style minimal model of kQ/(rels) by cells (s, t, W), W = path length.  Returns arr, d, wt, complete"""
    arr = {}; d = {}; wt = {}; LM = max(len(r) for r in rels)
    q = opp.Q(arr, [], 0); q.d = d
    cnt = 0; last_new = 1
    for W in range(1, WMAX + 1):
        if W == 1:
            for a, (s, t) in sorted(arrows.items()): arr[a] = (s, t, 0); wt[a] = 1
            continue
        for s in range(1, n + 1):
            bycell = {}
            def rec(v, w, cur):
                if w == W and len(cur) >= 2: bycell.setdefault(v, []).append(tuple(cur)); return
                for L, (a, b, _) in sorted(arr.items()):
                    if a == v and w + wt[L] <= W and wt[L] < W: rec(b, w + wt[L], cur + [L])
            rec(s, 0, [])
            for t in sorted(bycell):
                bydeg = {}
                for w in bycell[t]: bydeg.setdefault(q.deg(w), []).append(w)
                for m in sorted(bydeg):
                    Z = P.kernel(bydeg[m], [q.dpath(w) for w in bydeg[m]]) if m > 0 else [{w: 1} for w in bydeg[m]]
                    B = {}
                    for w in bydeg.get(m + 1, []): P.rref_add(B, q.dpath(w))
                    if m == 0:   # degree-0 words are arrow paths: the nonzero ones must survive as H_0 = A
                        for w in bydeg[0]:
                            if nz(w, rels): P.rref_add(B, {w: 1})
                    reps = [z for z in Z if P.rref_add(B, z)]
                    for z in reps:
                        cnt += 1; L = ('r%d' if m == 0 else 'c%d_%%d' % (m + 1)) % cnt
                        arr[L] = (s, t, m + 1); d[L] = P.integral(z); wt[L] = W; last_new = W
        if W >= last_new + LM or sum(1 for v in arr.values() if v[2]) > cap: break
    complete = W >= last_new + LM and sum(1 for v in arr.values() if v[2]) <= cap
    if plant:
        L = next(L for L in d if len(d[L]) >= 2); d[L].pop(next(iter(d[L])))
    return arr, d, wt, complete

def check_mm(n, rels, arr, d, wt, WMAX):
    q = opp.Q(arr, [], 0); q.d = d
    bad = [L for L in d if q.dlin(d[L])]
    for s in range(1, n + 1):
        cells = {}
        def rec(v, w, cur):
            if cur: cells.setdefault((v, w), []).append(tuple(cur))
            for L, (a, b, _) in arr.items():
                if a == v and w + wt[L] <= WMAX: rec(b, w + wt[L], cur + [L])
        rec(s, 0, [])
        for (t, W), ws in cells.items():
            bydeg = {}
            for w in ws: bydeg.setdefault(q.deg(w), []).append(w)
            for m in sorted(set(bydeg) | {0}):
                z = len(P.kernel(bydeg.get(m, []), [q.dpath(w) for w in bydeg.get(m, [])])) if m > 0 else len(bydeg.get(0, []))
                Bb = {}
                for w in bydeg.get(m + 1, []): P.rref_add(Bb, q.dpath(w))
                want = sum(1 for w in bydeg.get(0, []) if nz(w, rels)) if m == 0 else 0
                if z - len(Bb) != want: bad.append((s, t, W, m, z - len(Bb), want))
    return bad

def run_checks(n, arrows, rels, i, arr, d, wt, WENG):
    q, m, letters, wtf, dval, theta, A = CP.kM(arr, d, wt, i)
    C, sw = CP.shifted(arr, wt, i, m, letters, wtf)
    if min(sw.values()) <= 0: return dict(C=C, minshift='<= 0, engine not run')
    out = CP.analyse(n, rels, arr, d, wt, i, WENG)
    out['loop_letters_parent'] = sorted(L for L, (s, t, g) in arr.items() if s == t)
    return out

def daudit(n, ar, arr, d, wt, i, Dover=None):
    """review item 3: do (H1w) and (H3) see only letters with source in D(i)?  Checks (1) the D(i)-sourced letters of kM are
    closed under d; (2) every Theta letter, every linear-part letter of a Theta letter and every lead word is D(i)-sourced;
    (3) (H1w), Theta' and the lead words computed on the D(i)-sourced letters alone equal those on all letters.
    Returns (closed, inD, same, number of non-D letters)."""
    q, m, letters, wtf, dval, theta, A = CP.kM(arr, d, wt, i)
    D = Dover or reach(n, ar, i); inD = lambda L: m.src(L) in D
    LD = [L for L in letters if inD(L)]
    closed = all(inD(x) for L in LD for w in dval[L] for x in w)
    lin = {x for L in theta for x in P.h.lin_part(dval, L)}
    def run(Ls):
        dv = {L: dval[L] for L in Ls}; th = {L for L in theta if L in dv}
        bad, piv = P.h1w(Ls, dv, th)
        if bad: return ('H1w bad', sorted(map(h.nm, bad)))
        H, dH, pi_ok, Th, found, why = P.part_shape(m, Ls, dv, th, piv)
        return (sorted(map(h.nm, Th)), sorted((h.nm(t), tuple(map(h.nm, w))) for t, w in (found[1].items() if found else [])), bool(found) or not Th)
    full = run(letters)
    leadD = all(inD(x) for t, w in (P.part_shape(m, letters, dval, theta, P.h1w(letters, dval, theta)[1])[4] or [None, {}])[1].items() for x in w) \
        if not P.h1w(letters, dval, theta)[0] else True
    allD = all(inD(L) for L in theta) and all(inD(x) for x in lin) and leadD
    same = closed and run(LD) == full
    return closed, allD, same, len(letters) - len(LD)

# ---------------------------------------------------------------- enumeration
def canon(n, E, rels_v):
    """canonical form under vertex relabelling; E = set of (s, t), rels_v = set of vertex paths"""
    best = None
    for pm in itertools.permutations(range(1, n + 1)):
        f = lambda v: pm[v - 1]
        key = (tuple(sorted((f(s), f(t)) for s, t in E)), tuple(sorted(tuple(f(v) for v in r) for r in rels_v)))
        if best is None or key < best: best = key
    return best

def connected(n, E):
    adj = {v: set() for v in range(1, n + 1)}
    for s, t in E: adj[s].add(t); adj[t].add(s)
    seen, st = {1}, [1]
    while st:
        v = st.pop()
        for u in adj[v]:
            if u not in seen: seen.add(u); st.append(u)
    return len(seen) == n

def candidate_relations(arrows, LREL):
    """paths of length 2..LREL that are vertex paths of a simple digraph (no parallel arrows), any (cycles allowed)"""
    out = []
    for L in range(2, LREL + 1):
        def rec(p):
            if len(p) == L: out.append(tuple(p)); return
            t = arrows[p[-1]][1]
            for b, (s, _) in sorted(arrows.items()):
                if s == t: rec(p + [b])
        for a in sorted(arrows): rec([a])
    return out

def antichain(rels):
    return not any(r1 != r2 and not nz(r2, [r1]) for r1 in rels for r2 in rels)

def enumerate_parents(n, MAXA, LREL, MAXR=4):
    pairs = [(s, t) for s in range(1, n + 1) for t in range(1, n + 1) if s != t]
    seen = set()
    for k in range(n - 1, MAXA + 1):
        for E in itertools.combinations(pairs, k):
            if not connected(n, E): continue
            arrows = {'x%d' % (j + 1): e for j, e in enumerate(E)}
            if not has_cycle(set(range(1, n + 1)), arrows): continue
            # need some vertex i with an out-arrow and acyclic D(i)
            if not any(any(s == i for s, _ in E) and not has_cycle(reach(n, arrows, i), arrows) for i in range(1, n + 1)): continue
            cand = candidate_relations(arrows, LREL)
            for r in range(1, MAXR + 1):
                for rels in itertools.combinations(cand, r):
                    if not antichain(rels): continue
                    paths = nonzero_paths(arrows, rels, maxlen=4 * n)
                    if paths is None or not gldim_finite(arrows, rels): continue
                    vp = lambda path: tuple([arrows[path[0]][0]] + [arrows[a][1] for a in path])
                    key = canon(n, E, {vp(x) for x in rels})
                    if key in seen: continue
                    seen.add(key)
                    yield arrows, list(rels), paths

def instances(n, arrows, rels, paths):
    """(side, i, arrows', rels', paths', exact gate, library gate, library J, cycle reaches i)"""
    for side in ('A', 'Aop'):
        ar, rl = (arrows, rels) if side == 'A' else opposite(arrows, rels)
        ps = paths if side == 'A' else [tuple(reversed(p)) for p in paths]
        for i in range(1, n + 1):
            if not any(s == i for s, _ in ar.values()): continue
            if has_cycle(reach(n, ar, i), ar): continue
            g, why = exact_gate(ar, rl, ps, i)
            lg, lj = library(n, ar, rl, i)
            yield side, i, ar, rl, ps, g, why, lg, lj

def on_cycle(n, ar, v):
    return any(s == v and v in reach(n, ar, t) for s, t in ar.values())

def cycle_upstream(n, ar, i):
    """some vertex on an oriented cycle has a path to i"""
    return any(on_cycle(n, ar, v) and i in reach(n, ar, v) for v in range(1, n + 1))

def fmt(arrows, rels):
    return ' '.join(f'{a}:{s}>{t}' for a, (s, t) in sorted(arrows.items(), key=lambda x: (len(x[0]), x[0]))) + ' | rels ' + \
           ' '.join('.'.join(r) for r in rels)

def ctl_kupisch(N=5, LMAX=4):
    """[GLDIM control] the Anick tail-graph test against 126's Kupisch syzygy test on every antichain of zero relations on
    the cyclic quiver Z_n, n = 2..N, lengths 2..LMAX (both sides of the comparison say 'infinite' on many classes)"""
    agree = dis = inf = 0
    for nn in range(2, N + 1):
        arrows = {'x%d' % s0: (s0, s0 % nn + 1) for s0 in range(1, nn + 1)}
        for ls in itertools.product([0] + list(range(2, LMAX + 1)), repeat=nn):
            rl = [(s0 + 1, l) for s0, l in enumerate(ls) if l]
            if not rl: continue
            rels = [tuple('x%d' % ((s0 - 1 + k) % nn + 1) for k in range(l)) for s0, l in rl]
            if not antichain(rels): continue
            a, b = gldim_finite(arrows, rels), CP.gldim_finite(CP.kupisch(nn, rl))
            agree += a == b; dis += a != b; inf += not b
    print(f'[GLDIM control] Anick tail test vs Kupisch test on Z_n antichains n <= {N}, length <= {LMAX}: agree {agree}, disagree {dis} (must be 0); infinite by Kupisch {inf}')

def seed():
    """n = 4: x1: 1->2, x2: 2->1 (oriented 2-cycle), x3: 2->3, x4: 3->4; zero relations x1x2, x2x1; i = 3"""
    n = 4; arrows = {'x1': (1, 2), 'x2': (2, 1), 'x3': (2, 3), 'x4': (3, 4)}
    ctl_kupisch()
    for rels in ([('x1', 'x2')], [('x2', 'x1')], [('x1', 'x2'), ('x2', 'x1', 'x3')], [('x1', 'x2', 'x1')]):
        print(f'    gl.dim finite (Anick): {gldim_finite(arrows, rels)}')
        if not gldim_finite(arrows, rels): continue
        paths = nonzero_paths(arrows, rels)
        arr, d, wt, comp = minimal_model(n, arrows, rels)
        top = max(wt.values()); mm = check_mm(n, rels, arr, d, wt, top + 6)
        print(f'=== seed {fmt(arrows, rels)}: dim A = {len(paths) + n}, letters {len(arr)}, complete {comp}, max weight {top}, '
              f'[MM] H = A to W = {top + 6}: {"YES" if not mm else "NO " + str(mm[:3])}')
        print('    loop letters:', {L: (arr[L], h.show({tuple(w): c for w, c in d[L].items()})) for L in arr if arr[L][0] == arr[L][1]})
        for side, i, ar, rl, ps, g, why, lg, lj in instances(n, arrows, rels, paths):
            line = f'  {side} i={i} exact gate {int(g)} library gate {int(lg)} library J0 {int(lj)}'
            if g:
                a2, d2, w2, c2 = minimal_model(n, ar, rl) if side == 'Aop' else (arr, d, wt, comp)
                out = run_checks(n, ar, rl, i, a2, d2, w2, 2 * (max(w2.values()) + 6))
                line += ' ' + ' '.join(f'{k} {v}' for k, v in out.items())
            else:
                line += f' refused by path {why}'
            print(line)
    # controls
    rels = [('x1', 'x2'), ('x2', 'x1', 'x3')]
    arr, d, wt, _ = minimal_model(n, arrows, rels, plant=True)
    print(f'[MM control] seed with one d-term dropped: H = A failures {len(check_mm(n, rels, arr, d, wt, 10))} (must be > 0)')
    rels = [('x1', 'x2'), ('x3', 'x4')]
    g, why = exact_gate(arrows, rels, nonzero_paths(arrows, rels), 3)
    print(f'[GATE control] seed + relation x3x4: exact gate at 3 = {g} (must be False), library gate {library(n, arrows, rels, 3)[0]} (must be False)')
    # [BARE control]: the glued-vs-bare comparison must see a difference where there is one
    ga = {'x1': (1, 2), 'x2': (2, 3), 'x3': (3, 4), 'a': (5, 6), 'b': (6, 5), 'g': (6, 1)}
    gr = [('x1', 'x2', 'x3'), ('a', 'b')]
    A1, D1, W1, _ = minimal_model(6, ga, gr); o1 = canon_out(run_checks(6, ga, gr, 1, A1, D1, W1, 2 * (max(W1.values()) + 6)), A1, D1, W1)
    oc = CP.analyse(6, gr, A1, D1, W1, 1, 2 * (max(W1.values()) + 6), ctrl=True)
    print(f'[ENG / SH controls] glued A4 30 + 2-cycle, i = 1: ENG {oc["ENG"]}; closed ideal (stars) in place of I: {oc.get("ENGctrl_stars")} '
          f'(bad cells must be > 0); planted g: {oc.get("SHctrl")}')
    res = []
    for lr in ([('x1', 'x2', 'x3')], [('x1', 'x2')], [('x2', 'x3')]):
        la = {x: ga[x] for x in ('x1', 'x2', 'x3')}
        A2, D2, W2, _ = minimal_model(4, la, lr); res.append(o1 == canon_out(run_checks(4, la, lr, 1, A2, D2, W2, 2 * (max(W2.values()) + 6)), A2, D2, W2))
    print(f'[BARE control] glued A4 30 + 2-cycle at i = 1 {o1}: equal to bare A4 30 / A4 20 / A4 02: {res} (must be [True, False, False])')
    # [SIMPLE]: non-simple nonzero path x2.x1.x3 (2 -> 1 -> 2 -> 3) killed by the only out-arrow, simple paths all survive
    for rels in ([('x1', 'x2', 'x1'), ('x2', 'x1', 'x3', 'x4')], [('x1', 'x2'), ('x2', 'x1', 'x3', 'x4')]):
        ps = nonzero_paths(arrows, rels)
        g, why = exact_gate(arrows, rels, ps, 3); lg, lj = library(n, arrows, rels, 3)
        from quivermutation import procedure
        alg = path_algebra([list(arrows[x]) for x in sorted(arrows)], [[[arrows[r[0]][0]] + [arrows[x][1] for x in r]] for r in rels])
        try: procedure.mutateAtVertex(alg.quiver, procedure.relationsFrom(alg), 3); mut = 'runs (D(3) acyclic, not refused)'
        except Exception as e: mut = type(e).__name__
        print(f'[SIMPLE control] {fmt(arrows, rels)}: gl.dim finite {gldim_finite(arrows, rels)}, exact gate at 3 = {g} '
              f'(refusing path {why}), library gate {lg}, library J0 {lj}, mutateAtVertex(3): {mut}')

def sweep(n, MAXA, LREL, plan=False):
    t0 = time.time(); tally = dict(parents=0, inst=0, admitted=0, lib_disagree=0, lib_J_disagree=0, upstream=0, incomplete=0,
                                   run=0, allpass=0, loopletters=0, theta=0, thetaprime=0, thetaprime_up=0, J_disagree_on_admitted=0,
                                   admitted_Aop=0, up_Aop=0, aud_closed=0, aud_inD=0, aud_same=0, aud_nonD=0)
    fails = []; disagree = []; aop_ex = []
    for arrows, rels, paths in enumerate_parents(n, MAXA, LREL):
        tally['parents'] += 1
        if plan: continue
        models = {}
        for side, i, ar, rl, ps, g, why, lg, lj in instances(n, arrows, rels, paths):
            tally['inst'] += 1
            if lg != g: tally['lib_disagree'] += 1; disagree.append((fmt(ar, rl), side, i, g, lg, why))
            if lj != g: tally['lib_J_disagree'] += 1; tally['J_disagree_on_admitted'] += g
            if not g: continue
            tally['admitted'] += 1
            up = cycle_upstream(n, ar, i); tally['upstream'] += up
            tally['admitted_Aop'] += side == 'Aop'; tally['up_Aop'] += up and side == 'Aop'
            if side == 'Aop' and len(aop_ex) < 3: aop_ex.append((fmt(ar, rl), i, 'up', up))
            if side not in models: models[side] = minimal_model(n, ar, rl)
            arr, d, wt, comp = models[side]
            if not comp: tally['incomplete'] += 1; continue
            if 'mm' + side not in models:
                top = max(wt.values()); models['mm' + side] = check_mm(n, rl, arr, d, wt, top + 4)
            mm = models['mm' + side]
            out = run_checks(n, ar, rl, i, arr, d, wt, 2 * (max(wt.values()) + 6))
            tally['run'] += 1
            a1, a2, a3, a4 = daudit(n, ar, arr, d, wt, i)
            tally['aud_closed'] += a1; tally['aud_inD'] += a2; tally['aud_same'] += a3; tally['aud_nonD'] += a4 > 0
            tally['loopletters'] += bool(out.get('loop_letters_parent'))
            ok = (not mm and out.get('d2') == 0 and out.get('nonhom') == 0 and out.get('LF') and out.get('C')
                  and out.get('H1w') and out.get('PI', True) and out.get('SH', True) and out.get('ENG', (1,))[0] == 0)
            tally['allpass'] += bool(ok)
            if not ok: fails.append((fmt(ar, rl), side, i, up, len(mm), out))
            tally['theta'] += out.get('theta', 0) > 0
            tally['thetaprime'] += bool(out.get('thetaprime')); tally['thetaprime_up'] += bool(out.get('thetaprime')) and up
            if out.get('thetaprime'):
                print(f'  [Theta\'] {side} i={i} {fmt(ar, rl)} up {up} lead {out.get("lead")} ENG {out.get("ENG")} LF {out.get("LF")} cycles {out.get("cycles")}')
    print(f'[SWEEP n={n} arrows<={MAXA} rel length<={LREL}] {tally}  ({time.time() - t0:.0f} s)')
    for x in aop_ex: print('  [A^op admitted]', x)
    for x in disagree[:8]: print('  [LIB != EXACT gate]', x)
    for x in fails[:8]: print('  [FAIL]', x)

def bag():
    """review item 2: the referee's glued instance outside the glue list (relation b.a.g.x1 of length 4), and a control
    for the D(i) audit (D(i) replaced by {i} must break closure)"""
    arrows = {'x1': (1, 2), 'x2': (2, 3), 'a': (4, 5), 'b': (5, 4), 'g': (5, 1)}
    rels = [('a', 'b'), ('b', 'a', 'g', 'x1')]
    ps = nonzero_paths(arrows, rels); g, why = exact_gate(arrows, rels, ps, 1); lg, lj = library(5, arrows, rels, 1)
    print(f'[BAG] {fmt(arrows, rels)} at i = 1: gl.dim finite {gldim_finite(arrows, rels)}, D(1) acyclic {not has_cycle(reach(5, arrows, 1), arrows)}, '
          f'exact gate {g} (refusing path {why}), library gate {lg}, library J0 {lj}')
    arrows = {'x1': (1, 2), 'x2': (2, 1), 'x3': (2, 3), 'x4': (3, 4)}; rels = [('x1', 'x2')]
    arr, d, wt, _ = minimal_model(4, arrows, rels)
    print(f'[AUDIT control] seed x1.x2 = 0 at i = 3: audit with D(3) = {sorted(reach(4, arrows, 3))}: {daudit(4, arrows, arr, d, wt, 3)}; '
          f'with D replaced by {{1, 2, 3}} (not successor-closed): closed {daudit(4, arrows, arr, d, wt, 3, Dover={1, 2, 3})[0]} (must be False)')

def lna_rows(nL):
    """all antichains of zero relations on the line 1 -> ... -> nL (start s, length >= 2), as lists of (s, l)"""
    cand = [(s0, l) for s0 in range(1, nL) for l in range(2, nL - s0 + 1)]
    for r in range(0, nL):
        for rl in itertools.combinations(cand, r):
            if len({s0 for s0, _ in rl}) < r: continue
            if any(a != b and a[0] <= b[0] and b[0] + b[1] <= a[0] + a[1] for a in rl for b in rl): continue
            yield list(rl)

def glue_parents(NL, long=False):
    """an LNA on 1..nL (nL <= NL) with an oriented 2-cycle c1 <-> c2 (vertices nL+1, nL+2) glued upstream by g: c2 -> v;
    cycle relations one of ab, ba, aba, bab; at most one relation through g.  n = nL + 2 <= NL + 2."""
    for nL in range(2, NL + 1):
        for row in lna_rows(nL):
            for v in range(1, nL + 1):
                c1, c2 = nL + 1, nL + 2
                arrows = {'x%d' % j: (j, j + 1) for j in range(1, nL)}
                arrows.update({'a': (c1, c2), 'b': (c2, c1), 'g': (c2, v)})
                lrel = [tuple('x%d' % j for j in range(s0, s0 + l)) for s0, l in row]
                xv = ('x%d' % v,) if v < nL else ()
                xv2 = ('x%d' % v, 'x%d' % (v + 1)) if v + 1 < nL else ()
                glue = [(), ('a', 'g'), ('b', 'a', 'g')] + ([('g',) + xv, ('a', 'g') + xv] if xv else []) + ([('g',) + xv2] if xv2 else [])
                if long:   # review item 2: glue relations of length 4, the only ones that can kill a non-simple path into i
                    xv3 = ('x%d' % v, 'x%d' % (v + 1), 'x%d' % (v + 2)) if v + 2 < nL else ()
                    glue = [x for x in ([('b', 'a', 'g') + xv] if xv else []) + ([('a', 'g') + xv2] if xv2 else []) + ([('g',) + xv3] if xv3 else [])]
                for crel in ([('a', 'b')], [('b', 'a')], [('a', 'b', 'a')], [('b', 'a', 'b')]):
                    for gr in glue:
                        rels = lrel + crel + ([gr] if gr else [])
                        if not antichain(rels): continue
                        paths = nonzero_paths(arrows, rels, maxlen=4 * (nL + 2))
                        if paths is None or not gldim_finite(arrows, rels): continue
                        yield nL + 2, arrows, rels, paths

def canon_out(out, arr, d, wt):
    """Theta' and lead words with parent letter names r<k>, c<m>_<k> replaced by what they are (the relation path, or
    (source, target, degree, weight)), so two models with different letter numbering compare"""
    import re
    def tag(mo):
        L = mo.group(0)
        if L.startswith('r'): return 'R[' + '.'.join(next(iter(d[L]))) + ']'
        s0, t0, g0 = arr[L]; return f'C[{s0},{t0},{g0},{wt[L]}]'
    sub = lambda x: re.sub(r'\b(r\d+|c\d+_\d+)\b', tag, x)
    return (sorted(sub(t) for t in (out.get('thetaprime') or [])), sorted((sub(t), sub(w)) for t, w in (out.get('lead') or [])))

def glue(NL, plan=False, long=False, only=False):
    t0 = time.time(); tally = dict(parents=0, admitted=0, admitted_A=0, admitted_Aop=0, up_Aop=0, aud_closed=0, aud_inD=0, aud_same=0, aud_nonD=0, lib_gate_wrong=0, lib_J_wrong=0, incomplete=0, run=0, allpass=0,
                                   theta=0, thetaprime=0, loop_upstream=0, Jwrong_on_admitted=0,
                                   up=0, up_through=0, thetaprime_up=0, thetaprime_up_through=0, bare_same=0, bare_diff=0)
    fails = []; leads = {}; examples = []; bare = {}
    for n, arrows, rels, paths in glue_parents(NL, long):
        if only and n - 2 != NL: continue
        tally['parents'] += 1
        if plan: continue
        models = {}
        for side, i, ar, rl, ps, g, why, lg, lj in instances(n, arrows, rels, paths):
            if lg != g: tally['lib_gate_wrong'] += 1; examples.append(('gate', fmt(ar, rl), side, i, g, lg, why))
            if lj != g: tally['lib_J_wrong'] += 1; tally['Jwrong_on_admitted'] += g
            if not g: continue
            tally['admitted'] += 1; tally['admitted_' + side] += 1
            if side not in models:
                arr, d, wt, comp = minimal_model(n, ar, rl)
                models[side] = (arr, d, wt, comp, check_mm(n, rl, arr, d, wt, max(wt.values()) + 4) if comp else ['incomplete'])
            arr, d, wt, comp, mm = models[side]
            if not comp: tally['incomplete'] += 1; continue
            out = run_checks(n, ar, rl, i, arr, d, wt, 2 * (max(wt.values()) + 6))
            tally['run'] += 1
            a1, a2, a3, a4 = daudit(n, ar, arr, d, wt, i)
            tally['aud_closed'] += a1; tally['aud_inD'] += a2; tally['aud_same'] += a3; tally['aud_nonD'] += a4 > 0
            ok = (not mm and out.get('d2') == 0 and out.get('nonhom') == 0 and out.get('LF') and out.get('C')
                  and out.get('H1w') and out.get('PI', True) and out.get('SH', True) and out.get('ENG', (1,))[0] == 0)
            tally['allpass'] += bool(ok)
            tally['theta'] += out.get('theta', 0) > 0; tally['thetaprime'] += bool(out.get('thetaprime'))
            tally['loop_upstream'] += any(arr[L][0] == arr[L][1] and arr[L][2] >= 1 for L in arr)
            up = cycle_upstream(n, ar, i); thr = up and any(ar[x][1] == i for r in rl for x in r[:-1])
            tally['up'] += up; tally['up_through'] += thr; tally['up_Aop'] += up and side == 'Aop'
            if side == 'Aop' and len([e for e in examples if e[0] == 'aop']) < 3: examples.append(('aop', fmt(ar, rl), i, 'up', up, 'thetaprime', out.get('thetaprime')))
            tally['thetaprime_up'] += up and bool(out.get('thetaprime')); tally['thetaprime_up_through'] += thr and bool(out.get('thetaprime'))
            if not ok: fails.append((fmt(ar, rl), i, len(mm), out))
            if out.get('thetaprime'):
                for t, w in out['lead']: leads[w] = leads.get(w, 0) + 1
            if up and side == 'A':   # the same step on the bare LNA (cycle, glue arrow and glue relation removed): same Theta' and leads?
                nL = n - 2; la = {x: e for x, e in ar.items() if x[0] == 'x'}; lr = [r for r in rl if all(x[0] == 'x' for x in r)]
                key = (tuple(la.items()), tuple(lr), i)
                if key not in bare:
                    if not lr: bare[key] = ([], [])
                    else:
                        A2, D2, W2, _ = minimal_model(nL, la, lr)
                        o2 = run_checks(nL, la, lr, i, A2, D2, W2, 2 * (max(W2.values()) + 6))
                        bare[key] = canon_out(o2, A2, D2, W2)
                same = bare[key] == canon_out(out, arr, d, wt)
                tally['bare_same'] += same; tally['bare_diff'] += not same
                if not same: examples.append(('bare_diff', fmt(ar, rl), i, out.get('thetaprime'), out.get('lead'), bare[key]))
                if len(examples) < 400 and thr: examples.append(('through_i', fmt(ar, rl), i, out.get('thetaprime'), out.get('lead'), out['ENG'], 'LF', out['LF'], out.get('g=pathlength')))
    print(f'[GLUE nL{"=" if only else "<="}{NL}{" long glue relations" if long else ""}] {tally}  ({time.time() - t0:.0f} s)')
    print('  lead words (x\'.q or a*.p shapes):', dict(sorted(leads.items(), key=lambda x: -x[1])))
    for x in [e for e in examples if e[0] == 'gate'][:6]: print('  [LIB gate != EXACT]', x)
    for x in [e for e in examples if e[0] == 'through_i'][:12]: print('  [cycle reaches i, relation through i]', x)
    for x in [e for e in examples if e[0] == 'aop'][:3]: print('  [A^op admitted]', x)
    for x in [e for e in examples if e[0] == 'bare_diff'][:12]: print('  [BARE != GLUED]', x)
    for x in fails[:8]: print('  [FAIL]', x)

if __name__ == '__main__':
    if sys.argv[1] == 'seed': seed()
    elif sys.argv[1] == 'glue': glue(int(sys.argv[2]), plan='--plan' in sys.argv, long='--long' in sys.argv, only='--only' in sys.argv)
    elif sys.argv[1] == 'bag': bag()
    else: sweep(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), plan='--plan' in sys.argv)
