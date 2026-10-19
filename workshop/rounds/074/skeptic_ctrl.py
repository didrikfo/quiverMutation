"""Round 074 (skeptic), response to referee: (1) op/dual comparison of the 4 new iso classes with the 12 classes (and the 8 window ones);
(2) positive/negative controls of the presentation iso test.  Needs DIR from skeptic_failkids.py and an_summary.pkl from skeptic_failkids_an.py.
  skeptic_ctrl.py DIR"""
import sys, pickle, random, collections, itertools
sys.path.insert(0, '.')
import networkx as nx
from quivermutation import procedure, fingerprint, arrowPaths
D = sys.argv[1]
# presentation = (arrows list of (x,y,k), rels list of {path(tuple of arrows): coef})
def pres(a): return (sorted(a.quiver.edges(keys=True)), [dict(r) for r in procedure.relationsFrom(a)])
def struct(P):
    arrows, rels = P; g = nx.DiGraph()
    for v in {x for a in arrows for x in a[:2]}: g.add_node(('v', v), lab='v')
    for (x, y, k) in arrows:
        g.add_node(('a', x, y, k), lab='a'); g.add_edge(('a', x, y, k), ('v', x), lab='s'); g.add_edge(('a', x, y, k), ('v', y), lab='t')
    for i, rel in enumerate(rels):
        g.add_node(('r', i), lab='r')
        for j, (path, c) in enumerate(rel.items()):
            g.add_node(('p', i, j), lab='p%s_%d' % (abs(c), len(path))); g.add_edge(('r', i), ('p', i, j), lab='rp')
            for q, arr in enumerate(path): g.add_edge(('p', i, j), ('a',) + tuple(arr), lab='pos%d' % q)
    return g
nm = lambda a, b: a['lab'] == b['lab']
def iso(P, Q): return nx.is_isomorphic(struct(P), struct(Q), node_match=nm, edge_match=nm)
def signedIso(A, B, limit=20000):
    gA, gB = struct(A), struct(B)
    gm = nx.isomorphism.DiGraphMatcher(gB, gA, node_match=nm, edge_match=nm)
    arrowsA = A[0]
    class Q: pass
    target = fingerprint._relationsUnder(A[1], {}, arrowsA, gauge=True)
    for t, m in enumerate(gm.isomorphisms_iter()):
        if t >= limit: return None
        ren = {k[1:]: v[1:] for k, v in m.items() if k[0] == 'a'}
        if fingerprint._relationsUnder(B[1], ren, arrowsA, gauge=True) == target: return True
    return False
def op(P):
    arrows, rels = P
    return ([(y, x, k) for (x, y, k) in arrows],
            [{tuple((y, x, k) for (x, y, k) in reversed(path)): c for path, c in r.items()} for r in rels])
def relabel(P, rng):
    arrows, rels = P; V = sorted({x for a in arrows for x in a[:2]}); W = V[:]; rng.shuffle(W); vm = dict(zip(V, W))
    byb = collections.defaultdict(list)
    for (x, y, k) in arrows: byb[(x, y)].append(k)
    km = {}
    for (x, y), ks in byb.items(): p = ks[:]; rng.shuffle(p); km.update({(x, y, k): kk for k, kk in zip(ks, p)})
    am = {a: (vm[a[0]], vm[a[1]], km[a]) for a in arrows}
    sg = {a: rng.choice([1, -1]) for a in arrows}             # sign gauge: rescale arrow by +-1
    nr = []
    for r in rels:
        items = [(tuple(am[a] for a in path), c * __import__('math').prod(sg[a] for a in path)) for path, c in r.items()]; rng.shuffle(items); nr.append(dict(items))
    rng.shuffle(nr)
    return ([am[a] for a in arrows], nr)
# --- reload the 70 edges in the same order as skeptic_failkids_an.py
res = {}
for c in (1, 2):
    fails = []
    for s in range(4):
        t, f, stop = pickle.load(open('%s/c%d_s%d.pkl' % (D, c, s), 'rb')); fails += f
    seen = set(); edges = []
    for (d, pk, v, pa, ch) in fails:
        if pk is not None:
            if (pk, v) in seen: continue
            seen.add((pk, v))
        edges.append((d, pk, v, pa, ch))
    res[c] = edges
allE = [(c,) + e for c in (1, 2) for e in res[c]]
S = pickle.load(open(D + '/an_summary.pkl', 'rb'))
assert len(allE) == len(S['allE']) == 70 and [(e[0], e[1], e[3]) for e in allE] == S['allE'], 'order mismatch'
acl, onlyN, inE = S['acl'], S['onlyN'], S['inE151']
P = [pres(e[5]) for e in allE]
reps = [r[0] for r in acl]
isNew = [r in onlyN for r in acl]
win = [i for i, r in enumerate(acl) if not isNew[i]]; new = [i for i, r in enumerate(acl) if isNew[i]]
print('classes', len(acl), 'window', len(win), 'new', len(new))
# (1) op comparison: op of each new class rep vs every class rep (sign-blind iso; then signed iso when blind matches)
print('== op comparison: op(new rep) vs class reps')
nmatch = 0
for ni in new:
    O = op(P[reps[ni]]); row = []
    for ci in range(len(acl)):
        b = iso(O, P[reps[ci]])
        s = signedIso(P[reps[ci]], O) if b else None
        row.append((ci, 'window' if not isNew[ci] else 'new', b, s))
        if b: nmatch += 1
    print(' new class', ni, 'edges', len(acl[ni]), 'arrows/rels', len(P[reps[ni]][0]), len(P[reps[ni]][1]), '-> op matches (blind, signed):', [r for r in row if r[2]] or 'none')
print('total op matches of new reps against all 12 reps:', nmatch)
print('== op sanity: op(op(X)) == X, and op of each window rep vs each window rep (self-duality / window duality, informational)')
print(' op op ok:', all(iso(op(op(P[r])), P[r]) and signedIso(P[r], op(op(P[r]))) for r in reps))
print(' window rep pairs matched by op:', [(a, b) for a in win for b in win if iso(op(P[reps[a]]), P[reps[b]])])
# (2) controls
rng = random.Random(74)
print('== positive control: random vertex+arrow-key relabelling, sign gauge, relation/term shuffles of window children')
ok = tot = 0
for ci in win:
    for i in rng.sample(acl[ci], min(len(acl[ci]), 3)):
        for t in range(3):
            Q = relabel(P[i], rng); tot += 1
            # copy must also be recognised by the sign-aware test and have all canonicalKey-relevant structure; check blind and signed
            b = iso(P[reps[ci]], Q); s = signedIso(P[reps[ci]], Q) if b else False
            ok += bool(b and s)
print(' relabelled copies recognised in own class (blind and signed):', ok, 'of', tot)
print(' also each tested copy vs. all OTHER 11 class reps (should be 0 matches):', sum(iso(relabel(P[reps[ci]], rng), P[reps[cj]]) for ci in range(len(acl)) for cj in range(len(acl)) if ci != cj))
print('== negative control: sign change that is NOT a gauge')
# flip the sign of one term of a relation with >1 term in a window rep: blind iso must stay True, signed must become False when no gauge undoes it
cnt = 0; sep = 0
for ci in range(len(acl)):
    A = P[reps[ci]]
    for ri, r in enumerate(A[1]):
        if len(r) > 1:
            p0 = next(iter(r)); B = (A[0], [dict(x) for x in A[1]]); B[1][ri][p0] = -B[1][ri][p0]
            cnt += 1; sep += (not signedIso(A, B)) and iso(A, B); break
print(' reps with a multi-term relation, sign-flipped one term:', cnt, '; blind-iso True but signed-iso False:', sep)
print('== negative control: pairs of different classes with same (arrows, relations) counts separated by exact iso')
pairs = [(a, b) for a in range(len(acl)) for b in range(a + 1, len(acl)) if (len(P[reps[a]][0]), len(P[reps[a]][1])) == (len(P[reps[b]][0]), len(P[reps[b]][1]))]
print(' pairs', len(pairs), 'isomorphic among them:', sum(iso(P[reps[a]], P[reps[b]]) for a, b in pairs), '; of those pairs, same WL hash:',
      sum(nx.weisfeiler_lehman_graph_hash(struct(P[reps[a]]), node_attr='lab', edge_attr='lab') == nx.weisfeiler_lehman_graph_hash(struct(P[reps[b]]), node_attr='lab', edge_attr='lab') for a, b in pairs))
