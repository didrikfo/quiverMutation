"""Round 073 (toolsmith), T10/T5 item 1: an equivalence certificate with a NO branch for the endpoints (parent A, child B) of ONE key-kept J != 0 step.
  toolsmith_certificate.py in.pkl CLS I [MODE]     in.pkl = rounds/049/toolsmith_collect.py output (`7 2 20000 in.pkl 100`, ~280 s); run from the repo root.
Certificate for A ~ B (derived equivalent): a chain A -> ... -> M1 <- ... <- L -> ... -> M2 <- ... <- B of F/R tilting edges from the printed E-157 paths
(rounds/049/toolsmith_paths_logs.txt), L = a seed LNA/dual.  Each edge is accepted ONLY if
  (e1) the independent Hom test (rounds/050/skeptic_tilt.py: Hom(T,T[m]) = 0 for m = -1, +1; no tiltingPlus / perI / gate / key) passes,
  (e2) the edge's child is vertex-preserving and its Cartan matrix equals dim Hom(T_i,T_j) (End(T) has the child's dimension data),
and each meeting M1, M2 is accepted ONLY if the two algebras are isomorphic as quivers with relations by an explicit networkx isomorphism of the
structure graph (vertices, arrows, relation terms) -- NOT canonicalKey.  The verdict is YES (chain closes) or NO (some edge or meet fails).
Generation of K^b(proj) by T is NOT tested (E-168 proves it per loopless step).
MODE chain (default): replay + verdict.   MODE corrupt: every single-edge corruption of the chain (other vertex / other kind) and every swap of a meet
for a non-isomorphic same-key algebra; the verifier must say NO.   MODE edges: all (kind, v) steps from every chain node, Hom verdict against
ground truth that does not use the Hom test (key-moving child => certainly inequivalent => must reject; tiltingPlus-true key-kept => must accept)."""
import sys, re, pickle, time, copy
import numpy as np
import networkx as nx
A_ = sys.argv; PK, CLS, IDX = A_[1], int(A_[2]), int(A_[3]); MODE = A_[4] if len(A_) > 4 else 'chain'
sys.argv = ['x', PK, str(CLS), '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
from quivermutation import invariants

def struct(a):
    g = nx.DiGraph()
    for v in a.quiver.nodes(): g.add_node(('v', v), lab='v')
    for x, y in sorted(set(a.quiver.edges())): g.add_edge(('v', x), ('v', y), lab='a')
    # parallel arrows: record multiplicity on the edge label
    for x, y in set(a.quiver.edges()): g.edges[('v', x), ('v', y)]['lab'] = 'a%d' % a.quiver.number_of_edges(x, y)
    for i, rel in enumerate(procedure.relationsFrom(a)):
        g.add_node(('r', i), lab='r')
        for j, term in enumerate(rel):
            g.add_node(('p', i, j), lab='p%d' % len(term)); g.add_edge(('r', i), ('p', i, j), lab='rp')
            for k, arr in enumerate(term):
                g.add_edge(('p', i, j), ('v', arr[0]), lab='s%d' % k)
            g.add_edge(('p', i, j), ('v', term[-1][1]), lab='e')
    return g
def iso(x, y):
    """explicit isomorphism test of quivers with relations. CAVEAT: parallel arrows are only counted, and the arrow identity inside relation terms
    is collapsed to its endpoints, so for parallel-arrow bundles this is a necessary condition (reported separately)."""
    return nx.is_isomorphic(struct(x), struct(y), node_match=lambda a, b: a['lab'] == b['lab'], edge_match=lambda a, b: a['lab'] == b['lab'])
def parallel(x): return any(x.quiver.number_of_edges(a, b) > 1 for a, b in set(x.quiver.edges()))

def step(x, kind, v):
    a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
    V = sorted(x.quiver.nodes)
    if not mutation.mutationIsPossibleAtVertex(a, v): return None
    raw = mutation.quiverMutationAtVertex(a, v)
    c = reduction.reducePathAlgebra(raw)
    if sorted(c.vertices()) != V: return None
    y = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
    return y, a, c
def horn(x, kind, v):
    """independent verdict: (child, tiltOK, cartanOK, dimHom(T,T[-1]) total) or None if no vertex-preserving step."""
    r = step(x, kind, v)
    if r is None: return None
    y, a, c = r
    Vs, res = stepTest(a, v, -1)
    zero = not np.array(res[-1]).any() and not np.array(res[1]).any()
    C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
    cart = bool((np.array(res[0], dtype=int) == C.T).all())
    return y, zero, cart, int(np.abs(np.array(res[-1])).sum())

recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
rec = recs[IDX]; A0, B0, v0 = rec['parentObj'], rec['childObj'], rec['v']
seeds = classes[base]
paths = {}
SEC = None
for line in open('workshop/rounds/049/toolsmith_paths_logs.txt'):
    if line.startswith('=='): SEC = line.strip('= \n'); continue
    m = re.match(r'(paths|parents) (\d+) key (\w+) parentdepth (\d+) replay ok total (\d+) \| child->meet: (.*?) \| LNA/dual #(\d+) ->meet: (.*?) \d+s', line)
    if m and SEC.endswith('_c%d' % CLS) and int(m.group(2)) == IDX:
        paths[m.group(1)] = (m.group(6).split(), int(m.group(7)), m.group(8).split())
assert set(paths) == {'paths', 'parents'}, 'need both a parent path and a child path for this index: %s' % list(paths)
# NB the log's 'paths' = child side, 'parents' = parent side

from collections import Counter
CAUSE = Counter()
def run_side(x0, mv, seed, tamper=None):
    """returns (ok, nedges, hom_total_nonzero, detail). tamper = (side_edge_index, kind, v) replaces that edge's move."""
    ok = True; ne = 0; det = []
    x = x0
    for idx, s_ in enumerate(mv):
        kind, v = s_[0], int(s_[1:])
        if tamper and tamper[0] == idx: kind, v = tamper[1], tamper[2]
        h = horn(x, kind, v)
        if h is None: CAUSE['nostep'] += 1; return False, ne, ['NOSTEP'], None
        y, zero, cart, hm1 = h; ne += 1
        det.append('%s%d%s' % (kind, v, '' if zero and cart else '!'))
        if not (zero and cart): ok = False; CAUSE['hom' if not zero else 'cartan'] += 1
        x = y
    return ok, ne, det, x
def side(x0, mv_a, seed_idx, mv_b, tam_a=None, tam_b=None, seedObj=None):
    ok1, n1, d1, m1 = run_side(x0, mv_a, None, tam_a)
    L = seedObj if seedObj is not None else seeds[seed_idx]
    ok2, n2, d2, m2 = run_side(L, mv_b, None, tam_b)
    if m1 is None or m2 is None: return False, n1 + n2, d1 + ['||'] + d2, None, None
    meet = iso(m1, m2)
    if not meet: CAUSE['meet'] += 1
    return ok1 and ok2 and meet, n1 + n2, d1 + ['||'] + d2, meet, (parallel(m1) or parallel(m2))

def verdict(ptam=None, ctam=None, pseed=None, cseed=None):
    pm, ps, pl = paths['parents']; cm, cs, cl = paths['paths']
    r1 = side(A0, pm, ps, pl, ptam[0] if ptam else None, ptam[1] if ptam else None, pseed)
    r2 = side(B0, cm, cs, cl, ctam[0] if ctam else None, ctam[1] if ctam else None, cseed)
    return r1, r2

t0 = time.time()
pm, ps, pl = paths['parents']; cm, cs, cl = paths['paths']
print('record', IDX, 'class', CLS, 'step v =', v0, 'parent depth', rec['depth'], 'parent parallel arrows', parallel(A0), 'child parallel arrows', parallel(B0))
# the failing step itself
h = horn(A0, 'F', v0)
print('failing step F%d: child iso to recorded child %s; independent Hom test tilting=%s cartan=%s dim Hom(T,T[-1]) total=%d; recorded J=%s tiltingPlus=%s' %
      (v0, iso(h[0], B0), h[1], h[2], h[3], rec['J'], rec['tp']))
print('parent->LNA#%d: %s || %s ; child->LNA#%d: %s || %s' % (ps, ' '.join(pm), ' '.join(pl), cs, ' '.join(cm), ' '.join(cl)))
if MODE == 'chain':
    r1, r2 = verdict()
    for nm, r in (('parent side', r1), ('child side', r2)):
        print(nm, 'ok=%s edges=%d meet-isomorphic=%s parallel-at-meet=%s' % (r[0], r[1], r[3], r[4]), ' '.join(r[2]))
    print('PARENT ~ LNA#%d ~ CHILD certificate:' % ps, 'YES' if (r1[0] and r2[0] and ps == cs) else ('YES-via-two-LNAs(%d,%d)' % (ps, cs) if r1[0] and r2[0] else 'NO'),
          'total edges', r1[1] + r2[1], '%.0fs' % (time.time() - t0))
elif MODE == 'corrupt':
    base_ok = all(r[0] for r in verdict())
    print('uncorrupted chain accepted:', base_ok)
    V = sorted(A0.quiver.nodes)
    tot = acc = 0; accl = []
    for side_name, mv_a, mv_b in (('parent', pm, pl), ('child', cm, cl)):
        for which, mv in (('A', mv_a), ('B', mv_b)):
            for e in range(len(mv)):
                for kind in 'FR':
                    for v in V:
                        if (kind, v) == (mv[e][0], int(mv[e][1:])): continue
                        tam = ((e, kind, v), None) if which == 'A' else (None, (e, kind, v))
                        if side_name == 'parent': r = verdict(ptam=tam)[0]
                        else: r = verdict(ctam=tam)[1]
                        tot += 1
                        if r[0]: acc += 1; accl.append((side_name, which, e, kind, v))
    print('single-edge corruptions: %d, accepted (verifier says YES) %d' % (tot, acc), accl[:10], '%.0fs' % (time.time() - t0))
    print('rejection causes (a corruption can trip several; counts are per failing edge/meet):', dict(CAUSE))
    # swap the LNA end for each other seed of the class (with the genuine moves): meets must fail unless the seed is isomorphic
    sw = []
    for si, L in enumerate(seeds):
        r = verdict(pseed=L)[0]; sw.append((si, bool(r[0])))
    print('parent side with LNA#%d replaced by each seed #k (true #%d):' % (ps, ps), sw)
elif MODE == 'edges':
    # ground truth without the Hom test: key-moving child => inequivalent; tiltingPlus true + key kept => equivalent edge
    nodes = []; x = A0
    for kind_v, x0 in ((pm, A0), (pl, seeds[ps]), (cm, B0), (cl, seeds[cs])):
        x = x0; nodes.append(x)
        for s_ in kind_v:
            x = horn(x, s_[0], int(s_[1:]))[0]; nodes.append(x)
    cnt = dict(keymove_total=0, keymove_rejected=0, eq_total=0, eq_accepted=0, keykept_tpfalse=0, keykept_tpfalse_rejected=0, nonvertex=0)
    bad = []
    for x in nodes:
        for kind in 'FR':
            a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
            for v in sorted(x.quiver.nodes):
                if not mutation.mutationIsPossibleAtVertex(a, v): continue
                raw = mutation.quiverMutationAtVertex(a, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                hh = horn(x, kind, v)
                if hh is None: cnt['nonvertex'] += 1; continue
                y, zero, cart, hm1 = hh; acc_ = zero and cart
                c = y if kind == 'F' else pathAlgebra.dualPathAlgebra(y)
                kk = search._coxeterKeyOrNone(c) == base
                tp = bool(tiltingPlus(a.quiver, procedure.relationsFrom(a), v)) and not perI(a, v)
                if not kk: cnt['keymove_total'] += 1; cnt['keymove_rejected'] += (not acc_)
                elif tp: cnt['eq_total'] += 1; cnt['eq_accepted'] += acc_
                else: cnt['keykept_tpfalse'] += 1; cnt['keykept_tpfalse_rejected'] += (not acc_)
                if (not kk and acc_) or (kk and tp and not acc_): bad.append((kind, v))
    print('chain nodes', len(nodes), cnt, 'disagreements', bad, '%.0fs' % (time.time() - t0))
elif MODE == 'keymove':
    # ground truth independent of the Hom test: a step whose child has a different Coxeter key cannot be a derived equivalence.
    cnt = Counter(); wrong = []
    layer = list(seeds) + [r['parentObj'] for r in recs] + [r['childObj'] for r in recs]
    for depth in range(1):
        nxt = []
        for x in layer:
            for kind in 'FR':
                a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
                if not nx.is_directed_acyclic_graph(a.quiver): continue
                for v in sorted(x.quiver.nodes):
                    if not mutation.mutationIsPossibleAtVertex(a, v): continue
                    raw = mutation.quiverMutationAtVertex(a, v)
                    if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                    hh = horn(x, kind, v)
                    if hh is None: cnt['nonvertex'] += 1; continue
                    y, zero, cart, hm1 = hh
                    acc_ = zero and cart
                    kk = search._coxeterKeyOrNone(y) == base
                    cnt[('keykept' if kk else 'keymoving', 'accepted' if acc_ else 'rejected')] += 1
                    if not kk and acc_: wrong.append((depth, kind, v))
                    if kk and not acc_: cnt['keykept-rejected-Hm1=%d' % hm1] += 1
                    if kk and acc_ and depth == 0: nxt.append(y)
        layer = nxt[:60]
    print(dict(cnt), 'keymoving accepted (must be 0):', wrong, '%.0fs' % (time.time() - t0))
elif MODE == 'seedlink':
    # join two seed LNAs/duals a, b (argv 5, 6) of the class by F/R tilting moves; bucket by canonicalKey, then CONFIRM every edge with the Hom test and the meet by iso.
    a_i, b_i, D = int(A_[5]), int(A_[6]), int(A_[7])
    def bfs(root):
        seen = {fingerprint.canonicalKey(root): (root, None, None)}; fr = [root]; keys = {}
        for d in range(D):
            nx_ = []
            for x in fr:
                kx = fingerprint.canonicalKey(x)
                for kind, v, c in moves(x):
                    kc = fingerprint.canonicalKey(c)
                    if kc is None or kc in seen: continue
                    seen[kc] = (c, kx, (kind, v)); nx_.append(c)
            fr = nx_
        return seen
    Sa, Sb = bfs(seeds[a_i]), bfs(seeds[b_i])
    common = set(Sa) & set(Sb)
    print('seed', a_i, 'ball', len(Sa), 'seed', b_i, 'ball', len(Sb), 'common keys', len(common), '%.0fs' % (time.time() - t0))
    def trace(S, k):
        out = []
        while S[k][1] is not None: out.append(S[k][2]); k = S[k][1]
        return out[::-1]
    if common:
        k = sorted(common, key=lambda k: len(trace(Sa, k)) + len(trace(Sb, k)))[0]
        ma = ['%s%d' % m for m in trace(Sa, k)]; mb = ['%s%d' % m for m in trace(Sb, k)]
        ok1, n1, d1, m1 = run_side(seeds[a_i], ma, None); ok2, n2, d2, m2 = run_side(seeds[b_i], mb, None)
        print('path', ' '.join(ma), '||', ' '.join(mb), 'Hom-verified', ok1 and ok2, 'meet iso', iso(m1, m2), 'edges', n1 + n2, 'parallel', parallel(m1) or parallel(m2))
elif MODE == 'isoseeds':
    for i in range(len(seeds)):
        for j in range(i + 1, len(seeds)):
            if iso(seeds[i], seeds[j]): print('seeds', i, j, 'isomorphic as quivers with relations; canonicalKey equal:', fingerprint.canonicalKey(seeds[i]) == fingerprint.canonicalKey(seeds[j]))
elif MODE == 'seedlinkiso':
    # as seedlink, but the meet is found up to ISOMORPHISM (WL-hash bucket, then explicit iso), since canonicalKey is labelled (seeds i and its mirror
    # have different keys although isomorphic).  argv 5 6 7 = seed a, seed b, depth.
    a_i, b_i, D = int(A_[5]), int(A_[6]), int(A_[7])
    def sig(x):
        g = struct(x)
        return nx.weisfeiler_lehman_graph_hash(g, node_attr='lab', edge_attr='lab', iterations=3)
    def bfs(root):
        seen = {fingerprint.canonicalKey(root): (root, None, None, 0)}; fr = [root]
        for d in range(D):
            nx_ = []
            for x in fr:
                kx = fingerprint.canonicalKey(x)
                for kind, v, c in moves(x):
                    kc = fingerprint.canonicalKey(c)
                    if kc is None or kc in seen: continue
                    seen[kc] = (c, kx, (kind, v), d + 1); nx_.append(c)
            fr = nx_
        return seen
    def trace(S, k):
        out = []
        while S[k][1] is not None: out.append(S[k][2]); k = S[k][1]
        return out[::-1]
    Sa, Sb = bfs(seeds[a_i]), bfs(seeds[b_i])
    buck = {}
    for k, (x, _, _, d) in Sb.items(): buck.setdefault(sig(x), []).append(k)
    best = None
    for ka, (x, _, _, da) in Sa.items():
        for kb in buck.get(sig(x), []):
            if best is not None and da + Sb[kb][3] >= best[0]: continue
            if iso(x, Sb[kb][0]): best = (da + Sb[kb][3], ka, kb)
    print('seed', a_i, 'ball', len(Sa), 'seed', b_i, 'ball', len(Sb), 'iso-meet', best and best[0], '%.0fs' % (time.time() - t0))
    if best:
        ma = ['%s%d' % m for m in trace(Sa, best[1])]; mb = ['%s%d' % m for m in trace(Sb, best[2])]
        ok1, n1, d1, m1 = run_side(seeds[a_i], ma, None); ok2, n2, d2, m2 = run_side(seeds[b_i], mb, None)
        print('path', ' '.join(ma), '||', ' '.join(mb), 'Hom-verified', ok1 and ok2, 'meet iso', iso(m1, m2), 'edges', n1 + n2, 'parallel', parallel(m1) or parallel(m2))
elif MODE == 'f010':
    # NO control: the n = 9 F-010 pair (3060000 vs 3304000: same Coxeter key, certified inequivalent, E-057 via Ladkani Cor 3.15) and a YES control
    # (end of a random J = 0 walk of length 2D from 3060000).  Same meet-up-to-isomorphism search as seedlinkiso.  argv 5 = D, 6 = walk seed
    import random
    D = int(A_[5]); rnd = random.Random(int(A_[6]))
    a9 = nk.LinearNakayamaAlgebra.fromClassName('3060000'); b9 = nk.LinearNakayamaAlgebra.fromClassName('3304000')
    base = search._coxeterKeyOrNone(a9); assert base == search._coxeterKeyOrNone(b9)
    def sig(x): return nx.weisfeiler_lehman_graph_hash(struct(x), node_attr='lab', edge_attr='lab', iterations=3)
    def bfs(root):
        seen = {fingerprint.canonicalKey(root): (root, None, None, 0)}; fr = [root]
        for d in range(D):
            nx_ = []
            for x in fr:
                kx = fingerprint.canonicalKey(x)
                for kind, v, c in moves(x):
                    kc = fingerprint.canonicalKey(c)
                    if kc is None or kc in seen: continue
                    seen[kc] = (c, kx, (kind, v), d + 1); nx_.append(c)
            fr = nx_
        return seen
    def trace(S, k):
        out = []
        while S[k][1] is not None: out.append(S[k][2]); k = S[k][1]
        return out[::-1]
    def meet(x, y):
        Sa, Sb = bfs(x), bfs(y); buck = {}
        for k, (z, _, _, d) in Sb.items(): buck.setdefault(sig(z), []).append(k)
        best = None
        for ka, (z, _, _, da) in Sa.items():
            for kb in buck.get(sig(z), []):
                if best is not None and da + Sb[kb][3] >= best[0]: continue
                if iso(z, Sb[kb][0]): best = (da + Sb[kb][3], ka, kb)
        verified = None
        if best:
            ma = ['%s%d' % m for m in trace(Sa, best[1])]; mb = ['%s%d' % m for m in trace(Sb, best[2])]
            ok1, n1, d1, m1 = run_side(x, ma, None); ok2, n2, d2, m2 = run_side(y, mb, None)
            verified = (ok1 and ok2 and iso(m1, m2), n1 + n2)
        return len(Sa), len(Sb), best and best[0], verified
    # YES control: random walk of length 2D from a9
    x = a9; walk = []
    for _ in range(2 * D):
        mv = moves(x); k_, v_, c_ = rnd.choice(mv); walk.append('%s%d' % (k_, v_)); x = c_
    print('walk', ' '.join(walk))
    print('YES control (a9 vs end of its own walk, depth', D, '):', meet(a9, x), '%.0fs' % (time.time() - t0), flush=True)
    print('NO  control (F-010 pair 3060000 vs 3304000, depth', D, '):', meet(a9, b9), '%.0fs' % (time.time() - t0))
