"""Round 079 copy of round 077 skeptic_f010join.py (only change: walk:SEED:L:START). Round 077 (skeptic), T10: depth-bounded J = 0 join test of the key-equal n = 9 F-010 pair 3060000 vs 3304000 (inequivalent, Ladkani Cor 3.15, E-057).
Reuses the E-212 machinery (toolsmith_certificate.py f010 mode: F/R moves with gate + J = 0 + tiltingPlus + key kept, ball keyed by canonicalKey,
meet up to isomorphism by WL hash + networkx iso) but is standalone (no pickle) and RESUMABLE, so a 5 + 5 ball can be built in several < 10 min commands.
  skeptic_f010join.py side NAME ROOT DEPTH BUDGET_S CKPT   expand the ball of ROOT (class string, or 'walk:SEED:L' = end of an L-step random J=0 walk from 3060000)
                                                           until depth DEPTH done or BUDGET_S seconds spent; state saved to CKPT; prints size per depth.
  skeptic_f010join.py meet CKPT_A CKPT_B                   best iso meet between the COMPLETED levels of two balls, Hom-verified path replay if found.
  skeptic_f010join.py plan ROOT SECONDS                    expand ~SECONDS and print per-depth sizes and expansions/s (sizing only).
Run from the repository root."""
import sys, os, time, pickle, random
A_ = sys.argv; MODE = A_[1]
sys.path.insert(0, '.'); sys.argv = ['x', 'none', '2', '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
from quivermutation import invariants
import networkx as nx, numpy as np
t0 = time.time()
base = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3060000'))
assert base == search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3304000'))

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
def iso(x, y): return nx.is_isomorphic(struct(x), struct(y), node_match=lambda a, b: a['lab'] == b['lab'], edge_match=lambda a, b: a['lab'] == b['lab'])
def sig(x): return nx.weisfeiler_lehman_graph_hash(struct(x), node_attr='lab', edge_attr='lab', iterations=3)
def mkroot(name):
    if name.startswith('walk:'):
        _, sd, L, *rest = name.split(':'); rnd = random.Random(int(sd)); x = nk.LinearNakayamaAlgebra.fromClassName(rest[0] if rest else '3060000')   # walk:SEED:L[:START] (START default 3060000; round 079 added START)
        for _ in range(int(L)):
            mv = moves(x); k_, v_, c_ = rnd.choice(mv); x = c_
        return x
    return nk.LinearNakayamaAlgebra.fromClassName(name)

def expand(root, depth, budget, ck):
    """state: seen {key: (alg, parentkey, move, depth)}, frontier [keys at level d], nxt [keys level d+1], idx, d, nonekeys"""
    if ck and os.path.exists(ck): st = pickle.load(open(ck, 'rb'))
    else:
        k0 = fingerprint.canonicalKey(root); assert k0 is not None
        st = dict(seen={k0: (root, None, None, 0)}, frontier=[k0], nxt=[], idx=0, d=0, none=0, sizes=[1], secs=0.0, nexp=0)
    T = time.time(); s0 = dict(stats); st.setdefault('stats', dict(gate=0, J=0, tp=0, key=0, ok=0))
    while st['d'] < depth:
        while st['idx'] < len(st['frontier']):
            if time.time() - T > budget: break
            kx = st['frontier'][st['idx']]; x = st['seen'][kx][0]
            for kind, v, c in moves(x):
                kc = fingerprint.canonicalKey(c)
                if kc is None: st['none'] += 1; continue
                if kc in st['seen']: continue
                st['seen'][kc] = (c, kx, (kind, v), st['d'] + 1); st['nxt'].append(kc)
            st['idx'] += 1; st['nexp'] += 1
        if st['idx'] < len(st['frontier']): break
        st['d'] += 1; st['frontier'] = st['nxt']; st['nxt'] = []; st['idx'] = 0; st['sizes'].append(len(st['seen']))
    st['secs'] += time.time() - T
    for k_ in stats: st['stats'][k_] += stats[k_] - s0[k_]   # per-ball counts over all resumes: gate = candidate vertices, J = rejected because J != 0, tp = rejected by tiltingPlus, key = key moved, ok = kept
    if ck: pickle.dump(st, open(ck + '.tmp', 'wb')); os.replace(ck + '.tmp', ck)
    return st
def report(st):
    print('completed depth', st['d'], 'sizes(cumulative)', st['sizes'], 'in-progress level %d: %d/%d expanded, %d new' % (st['d'] + 1, st['idx'], len(st['frontier']), len(st['nxt'])),
          'expansions', st['nexp'], 'None-key children dropped', st['none'], 'secs %.0f' % st['secs'], 'move counts (gate/J!=0/tiltingPlus-fail/key-moved/kept)', st.get('stats', 'not recorded'), flush=True)
def trace(st, k):
    out = []
    while st['seen'][k][1] is not None: out.append(st['seen'][k][2]); k = st['seen'][k][1]
    return out[::-1]
def level(st, maxd): return {k: v for k, v in st['seen'].items() if v[3] <= maxd}

def replay(x0, mv):
    ok = True; x = x0
    for kind, v in mv:
        a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
        Vs, res = stepTest(a, v, -1)
        zero = not np.array(res[-1]).any() and not np.array(res[1]).any()
        raw = mutation.quiverMutationAtVertex(a, v); c = reduction.reducePathAlgebra(raw)
        C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
        ok = ok and zero and bool((np.array(res[0], dtype=int) == C.T).all())
        x = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
    return ok, x

if MODE == 'side':
    name, root, depth, budget, ck = A_[2], A_[3], int(A_[4]), float(A_[5]), A_[6]
    report(expand(mkroot(root), depth, budget, ck))
elif MODE == 'plan':
    st = expand(mkroot(A_[2]), 99, float(A_[3]), None); report(st)
    print('expansions/s %.2f' % (st['nexp'] / st['secs']))
elif MODE == 'meet':
    Sa, Sb = pickle.load(open(A_[2], 'rb')), pickle.load(open(A_[3], 'rb'))
    da, db = Sa['d'], Sb['d']; La, Lb = level(Sa, da), level(Sb, db)
    buck = {}
    for k, v in Lb.items(): buck.setdefault(sig(v[0]), []).append(k)
    best = None; ncmp = 0
    for ka, va in La.items():
        for kb in buck.get(sig(va[0]), []):
            if best is not None and va[3] + Lb[kb][3] >= best[0]: continue
            ncmp += 1
            if iso(va[0], Lb[kb][0]): best = (va[3] + Lb[kb][3], ka, kb)
    print('complete depths', da, db, 'ball sizes', len(La), len(Lb), 'WL-bucket iso tests', ncmp, 'meet total', best and best[0], '%.0fs' % (time.time() - t0))
    if best:
        ma, mb = trace(Sa, best[1]), trace(Sb, best[2])
        ok1, m1 = replay(Sa['seen'][next(iter(Sa['seen']))][0], ma); ok2, m2 = replay(Sb['seen'][next(iter(Sb['seen']))][0], mb)
        print('path', ma, '||', mb, 'Hom-verified', ok1 and ok2, 'meet iso', iso(m1, m2))
