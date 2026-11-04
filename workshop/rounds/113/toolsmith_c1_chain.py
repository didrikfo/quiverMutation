"""Round 113 (toolsmith), T10 / H-015: Hom-verified tilting chains parent ~ child for E-201's 13 key-kept J != 0, D = 0 steps (n = 7, key class 1).
  toolsmith_c1_chain.py list                  -- print the 13 steps from the walk checkpoint /tmp/skeptic_c1_13_ck1.pkl
  toolsmith_c1_chain.py chain I [DEPTH CAP]   -- bidirectional search for step I (0..12): parent A and child B joined by F/R vertex tilting edges
  toolsmith_c1_chain.py ctl I                 -- same chain, then planted corruptions (verifier must say NO)
Checkpoint: `timeout 10m .venv/bin/python -u workshop/rounds/069/skeptic_c1_13.py 1 540 23000` run twice (resumable; 2nd call prints the analysis).
SEARCH uses the repo's J = 0 / tiltingPlus / gate / key filters (only to find candidate edges).  VERIFICATION uses none of them: each edge is replayed
from scratch, accepted iff Hom_K(T,T[-1]) = Hom_K(T,T[1]) = 0 (rounds/050/skeptic_tilt.py, chain maps mod homotopy) and dim Hom(T_i,T_j) = Cartan of the
child; each end meeting must be a structure-graph isomorphism (not canonicalKey; not the Cartan matrix).  The step A -> B itself must FAIL the Hom test.
Generation of K^b(proj) by T is not tested (E-168 proves it per loopless step).  Run from the repository root."""
import sys, time, pickle, copy
import numpy as np, networkx as nx
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
from quivermutation import invariants
import sympy
MODE = ARGV[1]; n = 7; CLS = 1
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[CLS]

def struct(a):
    g = nx.DiGraph()
    for v in a.quiver.nodes(): g.add_node(('v', v), lab='v')
    for x, y in set(a.quiver.edges()): g.add_edge(('v', x), ('v', y), lab='a%d' % a.quiver.number_of_edges(x, y))
    for i, rel in enumerate(procedure.relationsFrom(a)):
        g.add_node(('r', i), lab='r')
        for j, term in enumerate(rel):
            g.add_node(('p', i, j), lab='p%d' % len(term)); g.add_edge(('r', i), ('p', i, j), lab='rp')
            for k, arr in enumerate(term): g.add_edge(('p', i, j), ('v', arr[0]), lab='s%d' % k)
            g.add_edge(('p', i, j), ('v', term[-1][1]), lab='e')
    return g
def iso(x, y):   # necessary test only when parallel arrows are present (arrow identity in relation terms is collapsed)
    return nx.is_isomorphic(struct(x), struct(y), node_match=lambda a, b: a['lab'] == b['lab'], edge_match=lambda a, b: a['lab'] == b['lab'])
def parallel(x): return any(x.quiver.number_of_edges(a, b) > 1 for a, b in set(x.quiver.edges()))
def dual(x): return pathAlgebra.dualPathAlgebra(x)

def step(x, kind, v):
    a = x if kind == 'F' else dual(x)
    V = sorted(x.quiver.nodes)
    if not mutation.mutationIsPossibleAtVertex(a, v): return None
    raw = mutation.quiverMutationAtVertex(a, v)
    c = reduction.reducePathAlgebra(raw)
    if sorted(c.vertices()) != V: return None
    return (c if kind == 'F' else dual(c)), a, c
def k0det(a, k, dup=None):
    """det of the K0 classes of the n summands of T in K0(proj A) = Z^V: [P_i] (i != k), [T_k] = [P_k] - sum_{k->h} [P_h]; dup=j replaces T_k by a second copy of P_j (control P2)."""
    V = sorted(a.quiver.nodes); ix = {u: i for i, u in enumerate(V)}; M = np.zeros((len(V), len(V)), dtype=object)
    for u in V:
        if u != k: M[ix[u], ix[u]] = 1
    if dup is not None: M[ix[k], ix[dup]] = 1
    else:
        M[ix[k], ix[k]] = 1
        for h in set(a.quiver.successors(k)): M[ix[k], ix[h]] -= a.quiver.number_of_edges(k, h)
    return int(sympy.Matrix(M.tolist()).det())
def horn(x, kind, v):
    """independent verdict (child, Hom(T,T[+-1]) = 0, Cartan(child) = dim Hom(T_i,T_j), total dim Hom(T,T[-1]))."""
    r = step(x, kind, v)
    if r is None: return None
    y, a, c = r
    Vs, res = stepTest(a, v, -1)
    zero = not np.array(res[-1]).any() and not np.array(res[1]).any()
    C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
    k0 = k0det(a, v) in (1, -1); loopfree = a.quiver.number_of_edges(v, v) == 0   # round 113 response: loop-free = E-168 generation condition
    return y, zero, bool((np.array(res[0], dtype=int) == C.T).all()), int(np.abs(np.array(res[-1])).sum()), k0, loopfree

def moves(alg):   # SEARCH filter only
    if list(nx.simple_cycles(alg.quiver)): return []
    V = sorted(alg.quiver.nodes); out = []
    for kind in 'FR':
        a = alg if kind == 'F' else dual(alg); rels = procedure.relationsFrom(a)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(a, v): continue
            if perI(a, v) or not tiltingPlus(a.quiver, rels, v): continue
            raw = mutation.quiverMutationAtVertex(a, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            c = reduction.reducePathAlgebra(raw)
            if sorted(c.vertices()) != V: continue
            if kind == 'R': c = dual(c)
            if search._coxeterKeyOrNone(c) != base: continue
            out.append((kind, v, c))
    return out

def detpoly(C):
    x = sympy.symbols('x'); return sympy.Poly((x * sympy.Matrix(C.tolist()) + sympy.Matrix(C.T.tolist())).det(), x)
def load_steps():
    seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
    Z = [r for r in byCart.values() if detpoly(r['CA']) == detpoly(r['CB'])]
    return nexp, Z
def key(x): return fingerprint.canonicalKey(x)

def search_chain(A, B, depth, cap, ck=None, dead=420):
    """bidirectional BFS (resumable: state pickled at ck when the `dead` second deadline passes); returns (movesA, movesB, meet-key) or None."""
    import os
    T0 = time.time()
    if ck and os.path.exists(ck): st = pickle.load(open(ck, 'rb')); print('resumed search', ck, flush=True)
    else: st = dict(S=[{key(A): (None, None, None, A)}, {key(B): (None, None, None, B)}], fr=[[A], [B]], nxt=[[], []], pos=[0, 0], d=0, s=0)
    S, fr = st['S'], st['fr']
    LIM = [int(x) for x in os.environ.get('LIM', '%d,%d' % (depth, depth)).split(',')]   # per-side depth limits (A,B); LIM=5,6 resumes a 5,5 checkpoint
    st.setdefault('lv', [st['d'] + (st['s'] == 1), st['d']])   # d = completed rounds in old checkpoints
    while st['lv'][0] < LIM[0] or st['lv'][1] < LIM[1]:
        s = st['s']
        if st['lv'][s] >= LIM[s]: st['s'] = 1 - s; continue
        while st['pos'][s] < len(fr[s]):
            if ck and time.time() - T0 > dead:
                pickle.dump(st, open(ck, 'wb')); return 'PAUSED', (len(S[0]), len(S[1]))
            alg = fr[s][st['pos'][s]]; st['pos'][s] += 1; ka = key(alg)
            for kind, v, c in moves(alg):
                kc = key(c)
                if kc is None or kc in S[s]: continue
                S[s][kc] = (ka, kind, v, c); st['nxt'][s].append(c)
                if kc in S[1 - s]: return trace(S, s, kc, A, B), (len(S[0]), len(S[1]))
        fr[s] = st['nxt'][s]; st['nxt'][s] = []; st['pos'][s] = 0
        if len(S[0]) + len(S[1]) > cap: return None, (len(S[0]), len(S[1]))
        st['lv'][s] += 1
        if s == 1: st['d'] += 1
        st['s'] = 1 - s
    return None, (len(S[0]), len(S[1]))
def trace(S, s, kc, A, B):
    def back(sd, k):
        mv = []
        while S[sd][k][0] is not None:
            pk, kind, v, c = S[sd][k]; mv.append((kind, v)); k = pk
        return mv[::-1]
    # move sequences from A (resp. B) to the meeting algebra; both sides store the meeting key (it exists in S[1-s] already)
    return back(0, kc), back(1, kc)

def replay(x0, mv, tamper=None):
    ok = True; x = x0; log = []; fails = set()
    for i, (kind, v) in enumerate(mv):
        if tamper and tamper[0] == i: kind, v = tamper[1], tamper[2]
        h = horn(x, kind, v)
        if h is None: fails.add('nostep'); return False, None, log + ['NOSTEP'], fails
        y, zero, cart, hm, k0, lf = h; log.append('%s%d%s' % (kind, v, '' if zero and cart and k0 and lf else '!')); ok &= (zero and cart and k0 and lf); x = y
        if not zero: fails.add('Hom')
        if not k0: fails.add('K0')
        if not cart: fails.add('Cartan')
        if not lf: fails.add('loop')
    return ok, x, log, fails
VERIFY_FAILS = set()
def verify(A, B, mvA, mvB, tA=None, tB=None, swapB=None):
    okA, MA, la, fa = replay(A, mvA, tA); okB, MB, lb, fb = replay(B if swapB is None else swapB, mvB, tB)
    meet = MA is not None and MB is not None and iso(MA, MB)
    VERIFY_FAILS.clear(); VERIFY_FAILS.update(fa | fb | (set() if meet else {'meet'}))
    return okA and okB and meet, la, lb, meet, (MA is not None and parallel(MA)) or (MB is not None and parallel(MB))

t0 = time.time()
nexp, Z = load_steps()
if MODE == 'list':
    print('checkpoint nexp', nexp, 'D=0 steps', len(Z), 'kept', sum(search._coxeterKeyOrNone(r['ch']) == base for r in Z))
    for i, r in enumerate(Z):
        print(i, 'v', r['v'], 'depth', r['depth'], 'nexp', r['nexp'], 'J', dict(r['J']), 'tp', r['tp'], 'parent parallel', parallel(r['alg']), 'child parallel', parallel(r['ch']), 'dim', int(r['CA'].sum()), int(r['CB'].sum()))
    sys.exit()
I = int(ARGV[2]); DEPTH = int(ARGV[3]) if len(ARGV) > 3 else 5; CAP = int(ARGV[4]) if len(ARGV) > 4 else 20000
r = Z[I]; A, B, v0 = r['alg'], r['ch'], r['v']
print('step', I, 'v', v0, 'parent parallel', parallel(A), 'child parallel', parallel(B), flush=True)
h = horn(A, 'F', v0)
print('failing step F%d: child iso to recorded child %s; Hom(T,T[+-1])=0 %s; Cartan=dimHom %s; dim Hom(T,T[-1]) = %d  (must be NOT tilting)' % (v0, iso(h[0], B), h[1], h[2], h[3]), flush=True)
if MODE == 'replay':   # replay I 0 0 'moveA...' 'moveB...' : verify given edge lists (no search)
    mk = lambda t: [(w[0], int(w[1:])) for w in t.split()]
    ok, la, lb, meet, par = verify(A, B, mk(ARGV[5]), mk(ARGV[6]))
    print('REPLAY step', I, 'A-side', ' '.join(la), '| B-side', ' '.join(lb), 'meet-isomorphic', meet, 'parallel-at-meet', par, 'VERDICT', 'YES' if ok else 'NO', 'failing', sorted(VERIFY_FAILS)); sys.exit()
res, sizes = search_chain(A, B, DEPTH, CAP, ck='/tmp/c1chain_ck%d.pkl' % I)
print('search sizes (A side, B side)', sizes, '%.0fs' % (time.time() - t0), flush=True)
if res == 'PAUSED': print('STEP', I, 'PAUSED (rerun the same command to resume)', sizes); sys.exit()
if res is None: print('STEP', I, 'UNDECIDED: no chain within depth', DEPTH, 'cap', CAP); sys.exit()
mvA, mvB = res
ok, la, lb, meet, par = verify(A, B, mvA, mvB)
print('A-side', ' '.join(la), '| B-side', ' '.join(lb), '| edges', len(mvA) + len(mvB), 'meet-isomorphic', meet, 'parallel-at-meet', par)
print('STEP', I, 'CHAIN VERDICT', 'YES' if ok else 'NO (verifier rejected)', '%.0fs' % (time.time() - t0))
if MODE == 'ctl':
    tot = acc = 0; accl = []; anyf = {}; solo = {}
    V = sorted(A.quiver.nodes)
    for side, mv in ((0, mvA), (1, mvB)):
        for e in range(len(mv)):
            for kind in 'FR':
                for v in V:
                    if (kind, v) == mv[e]: continue
                    tam = (e, kind, v)
                    o = verify(A, B, mvA, mvB, tam if side == 0 else None, tam if side == 1 else None)[0]
                    tot += 1; acc += o
                    if not o:
                        for t in VERIFY_FAILS: anyf[t] = anyf.get(t, 0) + 1
                        if len(VERIFY_FAILS) == 1: t, = VERIFY_FAILS; solo[t] = solo.get(t, 0) + 1
                    if o: accl.append((side, e, kind, v))
    print('single-edge corruptions %d, accepted %d' % (tot, acc), accl[:10])
    print('rejected (%d): failing tests (each corruption may fail several) %s; rejected by ONE test only %s' % (tot - acc, anyf, solo))
    # swap the end: replace B by a different child of the 13 (a different algebra) -- meet must fail
    sw = []; swf = {}
    for j, rj in enumerate(Z):
        if j == I: continue
        o = verify(A, B, mvA, mvB, swapB=rj['ch'])
        sw.append((j, bool(o[0]))); swf[tuple(sorted(VERIFY_FAILS))] = swf.get(tuple(sorted(VERIFY_FAILS)), 0) + 1
    print('B replaced by another step\'s child, verdict YES for:', [j for j, o in sw if o], 'of', len(sw), 'failing-test sets', swf, '%.0fs' % (time.time() - t0))
