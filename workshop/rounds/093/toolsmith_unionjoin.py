"""Round 093 (toolsmith): T10, join of 2400230 (F/R ball, radius 8) vs 3304000 (radius 6) on the UNION of the J = 0 key-kept F/R moves and the lnaMoves edges.
Balls are state pickles of workshop/rounds/087/toolsmith_d6scan.py `build` (DEPTH=8 for 2400230, DEPTH=6 for 3304000).
usage (repository root):
  toolsmith_unionjoin.py --plan
  toolsmith_unionjoin.py lna A.pkl B.pkl      LNA census of both balls; lnaMoves components touched; component intersection (the lnaMoves jump), per-depth
  toolsmith_unionjoin.py extra B.pkl          replay the children of E-244 indices 8427 and 8063 from 3304000 (iso test, LNA test, radius-1 J = 0 ball)
  toolsmith_unionjoin.py meet A.pkl B.pkl [--extra]   F/R-only iso meet of the two complete balls (E-226 form, baseline); --extra adds the 8427 child as a B node at depth 7
The lnaMoves jump is an over-approximation: an LNA node of a ball is allowed to jump to ANY member of its lnaMoves component at cost 1 (real cost >= 1)."""
import sys, time, pickle, collections
MYA = list(sys.argv); MODE = MYA[1] if len(MYA) > 1 else '--plan'
if MODE == '--plan':
    print("""sizing (round 093): A = 2400230 depth 8 (20 210 nodes, build ~5 min in two commands, E-254), B = 3304000 depth 6 (15 401 nodes, build ~11 min in two commands, E-244).
 lna: standardise 35 600 nodes (~ms each) + closureUnderMoves per distinct class name (<= 380 components at n = 9, E-238): expected < 3 min.
 extra: replay 2 paths of 7 steps + radius-1 ball: < 1 min.   meet: WL signature of 35 600 graphs + iso on buckets: ~ 5-8 min (E-254 meets ran in the same form).
 each command < 10 min; builds: DEPTH=8 ... build 2400230 A8.pkl --budget-hours 0.14 ; DEPTH=6 ... build 3304000 B6.pkl --budget-hours 0.14 (exit 2 = rerun)""")
    sys.exit(0)
sys.argv = ['x', 'ctl', 'x']
src = open('workshop/rounds/087/toolsmith_d6scan.py').read().split("if MODE == 'build'")[0]
sys.path.insert(0, '.'); exec(compile(src, 's', 'exec'))
MODE = MYA[1]
from quivermutation import lnaMoves, lines
_sj = open('workshop/rounds/082/skeptic_join.py').read(); exec(_sj[_sj.index('def child('):_sj.index('def mkroot(')])  # child, struct, iso, sig
import networkx as nx, numpy as np
T0 = time.time()
def path(S, k):
    p = []
    while S['seen'][k][1] is not None: p.append(S['seen'][k][2]); k = S['seen'][k][1]
    return p[::-1]
def lnaName(x):
    q = x.quiver
    if q.number_of_edges() != len(q.nodes) - 1 or any(q.out_degree(v) > 1 or q.in_degree(v) > 1 for v in q.nodes): return None
    try: _, rl, _ = lnaMoves.standardise(x)
    except Exception: return None
    return None if rl is None else lines.className(rl)
_comp = {}
def component(name):
    if name not in _comp:
        orb = lnaMoves.closureUnderMoves(9, [int(c) for c in name])      # forward closure from this row
        _comp[name] = frozenset(orb)
    return _comp[name]
def ucomp(names):
    """union-find over the FULL undirected lnaMoves graph: every row of length 9 (nk.allRelationLengths), edges = movesByRule names (as in round 085 skeptic_hom.py comp)"""
    par = {}
    def f(a):
        while par.setdefault(a, a) != a: par[a] = par[par[a]]; a = par[a]
        return a
    for r in nk.allRelationLengths(9):
        u = ''.join(map(str, r))
        for v in lnaMoves.movesByRule(9, [int(c) for c in u]): par[f(u)] = f(v)
    return f
def replayPath(root, mv):
    x = root
    for kind, v in mv: x = child(x, kind, v)
    return x
if MODE == 'lna':
    SA, SB = pickle.load(open(MYA[2], 'rb')), pickle.load(open(MYA[3], 'rb'))
    print('depths', SA['d'], SB['d'], 'sizes', len(SA['seen']), len(SB['seen']), flush=True)
    info = {}
    for tag, S in (('A', SA), ('B', SB)):
        c = collections.Counter(); names = collections.defaultdict(list)
        for k, (x, pk, mv, d) in S['seen'].items():
            nm = lnaName(x)
            if nm is not None: names[nm].append(d); c[d] += 1
        info[tag] = names
        print(tag, 'LNA nodes by depth', dict(sorted(c.items())), 'distinct LNA classes', len(names), flush=True)
        print(tag, 'LNA classes (name: min depth)', {n_: min(v) for n_, v in sorted(names.items(), key=lambda t: min(t[1]))}, flush=True)
    allnames = set(info['A']) | set(info['B'])
    f = ucomp(sorted(allnames))
    cA = collections.defaultdict(list); cB = collections.defaultdict(list)
    for nm, ds in info['A'].items(): cA[f(nm)].append((min(ds), nm))
    for nm, ds in info['B'].items(): cB[f(nm)].append((min(ds), nm))
    print('undirected lnaMoves components touched: A %d, B %d, common %d (time %.0fs)' % (len(cA), len(cB), len(set(cA) & set(cB)), time.time() - T0))
    sz = collections.Counter(f(''.join(map(str, r))) for r in nk.allRelationLengths(9))
    print('touched components, full size (rows): A', {r: (sz[r], sorted(cA[r])[:3]) for r in cA}, 'B', {r: (sz[r], sorted(cB[r])[:3]) for r in cB})
    for r in set(cA) & set(cB):
        print('COMMON component', r, 'A members (depth,name)', sorted(cA[r]), 'B members', sorted(cB[r]), 'sizes of forward closures', 'n/a')
    print('B component of 3304000:', sorted(m for m in allnames if f(m) == f('3304000')), 'A component of 2400230:', sorted(m for m in allnames if f(m) == f('2400230')))
    print('same component 2400230 / 3304000 in the lnaMoves graph restricted to what the balls touch:', f('2400230') == f('3304000'))
elif MODE == 'extra':
    SB = pickle.load(open(MYA[2], 'rb')); root = SB['seen'][next(iter(SB['seen']))][0]
    PTH = dict(i8427=[('F', 2), ('R', 9), ('R', 8), ('R', 7), ('R', 9), ('R', 7)], i8063=[('F', 2), ('R', 8), ('R', 7), ('R', 9), ('F', 7), ('R', 9)])
    St = dict(i8427=('R', 7), i8063=('R', 9)); C = {}
    for nm, mv in PTH.items():
        par = replayPath(root, mv); jc = jcands(par); print(nm, 'parent key-kept J != 0 candidates', [t for t in jc if t[2] == 'kept'], 'parent in ball B', fingerprint.canonicalKey(par) in SB['seen'])
        C[nm] = child(par, *St[nm]); print(nm, 'child LNA name', lnaName(C[nm]), 'key kept', search._coxeterKeyOrNone(C[nm]) == base, flush=True)
    print('children isomorphic (labelled quiver + relations, up to vertex relabelling):', iso(C['i8427'], C['i8063']))
    print('child 8427 labelled-equal to child 8063:', fingerprint.canonicalKey(C['i8427']) == fingerprint.canonicalKey(C['i8063']))
    nb = [(k, v, fingerprint.canonicalKey(c)) for k, v, c in moves(C['i8427'])]
    print('radius-1 J = 0 ball of child 8427:', len(nb), 'neighbours; any neighbour key in ball B:', sum(1 for _, _, kk in nb if kk in SB['seen']), 'any iso to a B level-6 node:', 'see meet --extra')
    pickle.dump(dict(c8427=C['i8427'], c8063=C['i8063']), open(MYA[3], 'wb')) if len(MYA) > 3 else None
elif MODE == 'meet':
    SA, SB = pickle.load(open(MYA[2], 'rb')), pickle.load(open(MYA[3], 'rb'))
    La = dict(SA['seen']); Lb = dict(SB['seen'])
    extra = {}
    if '--extra' in MYA:
        root = SB['seen'][next(iter(SB['seen']))][0]
        par = replayPath(root, [('F', 2), ('R', 9), ('R', 8), ('R', 7), ('R', 9), ('R', 7)]); c = child(par, 'R', 7)
        Lb[('extra', 0)] = (c, None, None, 7)
        for k, v, c2 in moves(c): Lb[('extra', 1, k, v)] = (c2, None, None, 8)
        print('B gains child of 8427 (depth 7) and its radius-1 ball (depth 8)', flush=True)
    buck = {}
    for k, v in Lb.items(): buck.setdefault(sig(v[0]), []).append(k)
    best = None; ncmp = 0; nb = 0
    for ka, va in La.items():
        for kb in buck.get(sig(va[0]), []):
            nb += 1
            if best is not None and va[3] + Lb[kb][3] >= best[0]: continue
            ncmp += 1
            if iso(va[0], Lb[kb][0]): best = (va[3] + Lb[kb][3], ka, kb)
    print('ball sizes', len(La), len(Lb), 'WL-bucket pairs', nb, 'iso tests', ncmp, 'meet total', best and best[0], 'time %.0fs' % (time.time() - T0))
