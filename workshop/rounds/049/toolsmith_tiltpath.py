"""Round 049 (toolsmith), T10 (i): tilting-only meet-in-the-middle path search from the E-152 children back to an LNA.
Moves (all J = 0, tiltingPlus true, gate-admitted, legal, vertex-preserving, Coxeter key kept):
  F v : the forward step at v;   R v : the inverse of a forward step, taken as the forward step at v of the OPPOSITE algebra, carried back.
Mixed F/R paths allowed (an undirected graph on algebras, keyed by fingerprint.canonicalKey; a node with no canonical key is not deduplicated).
Target ball: all LNAs and duals of the child's key class, expanded dT levels with the same moves.  Child ball: dC levels.  A shared key = a
tilting-only path child ~ LNA of length <= dC + dT (hit proves derived equivalence to that LNA, given F/R are tilting steps); a miss proves nothing.
Controls: random forward J = 0 walks of length L from class LNAs, searched identically (known path of length L exists).
Usage: toolsmith_tiltpath.py in.pkl cls dT dC capNodes mode   mode = children | ctrlN:L  (N controls of walk length L)   run from repo root."""
import sys, pickle, time, copy, random
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
pk, cls, dT, dC, cap, mode = ARGV[1], int(ARGV[2]), int(ARGV[3]), int(ARGV[4]), int(ARGV[5]), ARGV[6]
N = 7
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)): classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
stats = dict(gate=0, J=0, tp=0, key=0, ok=0)

def moves(alg):
    """children of alg under F and R, J = 0 and tiltingPlus, key kept."""
    if list(nx.simple_cycles(alg.quiver)): return []
    V = sorted(alg.quiver.nodes); out = []
    for kind in 'FR':
        a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
        rels = procedure.relationsFrom(a)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(a, v): continue
            stats['gate'] += 1
            if perI(a, v): stats['J'] += 1; continue
            if not tiltingPlus(a.quiver, rels, v): stats['tp'] += 1; continue
            raw = mutation.quiverMutationAtVertex(a, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            c = reduction.reducePathAlgebra(raw)
            if sorted(c.vertices()) != V: continue
            if kind == 'R': c = pathAlgebra.dualPathAlgebra(c)
            if search._coxeterKeyOrNone(c) != base: stats['key'] += 1; continue
            stats['ok'] += 1; out.append((kind, v, c))
    return out

def ball(starts, depth, capn, label):
    seen = {}; frontier = []
    for s in starts:
        k = fingerprint.canonicalKey(s)
        if k is not None and k in seen: continue
        if k is not None: seen[k] = 0
        frontier.append(s)
    levels = [len(frontier)]; capped = False
    for d in range(1, depth + 1):
        nxt = []
        for alg in frontier:
            if len(seen) >= capn: capped = True; break
            for kind, v, c in moves(alg):
                k = fingerprint.canonicalKey(c)
                if k is None: nxt.append(c); continue
                if k in seen: continue
                seen[k] = d; nxt.append(c)
        levels.append(len(nxt)); frontier = nxt
        if capped: break
    return seen, levels, capped

def reach(start, tball, depth, capn):
    """BFS from start; first level d at which a key lies in tball -> (d + tball depth, d, td), levels, capped."""
    k0 = fingerprint.canonicalKey(start)
    if k0 in tball: return (tball[k0], 0, tball[k0]), [1], False
    seen = {k0}; frontier = [start]; levels = [1]; best = None; capped = False
    for d in range(1, depth + 1):
        nxt = []
        for alg in frontier:
            if len(seen) >= capn: capped = True; break
            for kind, v, c in moves(alg):
                k = fingerprint.canonicalKey(c)
                if k is None: nxt.append(c); continue
                if k in seen: continue
                seen.add(k); nxt.append(c)
                if k in tball and (best is None or d + tball[k] < best[0]): best = (d + tball[k], d, tball[k])
        levels.append(len(nxt)); frontier = nxt
        if best is not None or capped: break   # first level with a hit; total = d + tball depth (min over this level)
    return best, levels, capped

t0 = time.time()
import os, hashlib
cache = '/tmp/tsm/tball_c%d_d%d.pkl' % (cls, dT)   # cache of the target ball (scratch only; delete to recompute)
if os.path.exists(cache): tball, tl, tc = pickle.load(open(cache, 'rb'))
else:
    tball, tl, tc = ball(classes[base], dT, cap, 'target'); pickle.dump((tball, tl, tc), open(cache, 'wb'))
def tag(alg): return hashlib.md5(repr(fingerprint.canonicalKey(alg)).encode()).hexdigest()[:6]
print('target ball: class', cls, 'seeds', len(classes[base]), 'levels', tl, 'keys', len(tball), 'capped', tc, '%.0fs' % (time.time() - t0), flush=True)

def run(name, start, extra=''):
    t1 = time.time(); best, lv, capd = reach(start, tball, dC, cap)
    print(name, 'key', tag(start), extra, 'HIT total %d (child side %d + target side %d)' % best if best else 'MISS', 'child-ball levels', lv, 'capped' if capd else 'uncapped', '%.0fs' % (time.time() - t1), flush=True)
    return best

if mode == 'children':
    recs = [r for r in pickle.load(open(pk, 'rb')) if r['kind'] == 'fail']
    print('children', len(recs), flush=True); hits = 0
    for i, r in enumerate(recs):
        if run('child %d' % i, r['childObj'], 'parentdepth %d v %d' % (r['depth'], r['v'])): hits += 1
    print('SUMMARY cls', cls, 'children', len(recs), 'hits', hits, 'dT', dT, 'dC', dC, 'cap', cap, 'stats', stats)
else:
    nc, L = [int(x) for x in mode[4:].split(':')]; rnd = random.Random(1234 + L); hits = 0
    for i in range(nc):
        alg = copy.deepcopy(rnd.choice(classes[base])); path = []
        for step in range(L):
            ms = [m for m in moves(alg) if m[0] == 'F']
            if not ms: break
            kind, v, alg = rnd.choice(ms); path.append(v)
        if len(path) < L: print('control', i, 'walk stuck at', len(path)); continue
        if run('control %d' % i, alg, 'forward J=0 walk length %d path %s' % (L, path)): hits += 1
    print('SUMMARY controls L', L, 'n', nc, 'hits', hits, 'dT', dT, 'dC', dC, 'cap', cap, 'stats', stats)
