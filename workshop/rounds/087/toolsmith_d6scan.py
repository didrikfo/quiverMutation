"""Round 087 (toolsmith): T10, depth-6 scan of the 3304000 ball for key-kept perI (J != 0) steps. Checkpointed.
usage (repository root; env KEYROOT default 3060000 = same Coxeter key as 3304000, DEPTH default 6):
  toolsmith_d6scan.py --plan                          print the sizing model (no computation)
  toolsmith_d6scan.py build ROOT STATE.pkl [--budget-hours H]   resumable BFS; checkpoint every 300 expansions; exit 2 if budget spent (rerun same command to resume), 0 when done
  toolsmith_d6scan.py parents STATE.pkl PAR.pkl       write the depth-DEPTH parents (node, path) after build finished
  toolsmith_d6scan.py scan PAR.pkl START END [--budget-hours H]  scan parents[START:END]; prints KEPT lines, 'DONE upto i' ; exit 2 if budget spent (reports resume index)
"""
import sys, time, pickle, collections, os
A_ = [a for a in sys.argv]; MODE = A_[1] if len(A_) > 1 else '--plan'
BUD = 1e9
if '--budget-hours' in A_:
    i = A_.index('--budget-hours'); BUD = float(A_[i + 1]) * 3600; del A_[i:i + 2]
T00 = time.time()
if MODE == '--plan':
    print("""sizing model (3304000, from E-239): ball sizes 1,11,64,282,1088,4046 (depth 5); expansion ~0.08 s/node (116 s for the 1446 nodes of depth 0-4);
 depth-6 new nodes ~ 3.7 x 2958 = ~11 000; build of depth 6 = expand 2 958 depth-5 nodes ~ 240 s (+116 s if from scratch) => ~6 min total
 scan ~35 s per 1500 parents => ~11 000 parents ~ 4.5 min (shards of 4000: ~2 min each)
 3060000: depth-5 front 4 770 -> expand ~380 s, depth-6 parents ~ 17 000-19 000, scan ~7 min; build from scratch ~10 min -> use checkpointed build (2 commands)
 each command: add --budget-hours 0.14 (8.4 min); exit 2 = rerun same command""")
    sys.exit(0)
sys.path.insert(0, '.'); sys.argv = ['x', 'none', '2', '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
KEYROOT = os.environ.get('KEYROOT', '3060000'); DEPTH = int(os.environ.get('DEPTH', '6'))
base = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName(KEYROOT))
def jcands(alg):
    V = sorted(alg.quiver.nodes); out = []
    for kind in 'FR':
        a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(a, v) or not perI(a, v): continue
            try:
                raw = mutation.quiverMutationAtVertex(a, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): out.append((kind, v, 'illegal')); continue
                c = reduction.reducePathAlgebra(raw)
                if sorted(c.vertices()) != V: out.append((kind, v, 'verts')); continue
                if kind == 'R': c = pathAlgebra.dualPathAlgebra(c)
                k = search._coxeterKeyOrNone(c)
                out.append((kind, v, 'none' if k is None else ('kept' if k == base else 'moved')))
            except Exception as e: out.append((kind, v, 'err ' + type(e).__name__))
    return out
def rels(a): return [list(r) for r in procedure.relationsFrom(a)]
def left(): return BUD - (time.time() - T00)
if MODE == 'build':
    ROOT, OUT = A_[2], A_[3]
    if os.path.exists(OUT):
        S = pickle.load(open(OUT, 'rb')); print('resume depth', S['d'], 'idx', S['i'], 'ball', len(S['seen']), flush=True)
    else:
        root = nk.LinearNakayamaAlgebra.fromClassName(ROOT)
        S = dict(seen={fingerprint.canonicalKey(root): (root, None, None, 0)}, front=[root], nxt=[], d=0, i=0)
    def save():
        pickle.dump(S, open(OUT + '.tmp', 'wb')); os.replace(OUT + '.tmp', OUT)
    since = 0
    while S['d'] < DEPTH:
        front = S['front']
        while S['i'] < len(front):
            x = front[S['i']]; kx = fingerprint.canonicalKey(x)
            for kind, v, c in moves(x):
                k = fingerprint.canonicalKey(c)
                if k is not None and k not in S['seen']: S['seen'][k] = (c, kx, (kind, v), S['d'] + 1); S['nxt'].append(c)
            S['i'] += 1; since += 1
            if since >= 300 or left() < 20:
                save(); since = 0
                if left() < 20: print('BUDGET spent at depth', S['d'], 'idx', S['i'], '/', len(front), 'ball', len(S['seen'])); sys.exit(2)
        S['d'] += 1; S['front'] = S['nxt']; S['nxt'] = []; S['i'] = 0
        print('depth', S['d'], 'ball', len(S['seen']), '%.0fs' % (time.time() - T00), flush=True); save()
    print('BUILD DONE depth', DEPTH, 'ball', len(S['seen'])); sys.exit(0)
elif MODE == 'parents':
    S = pickle.load(open(A_[2], 'rb')); seen = S['seen']
    assert S['d'] == DEPTH
    def path(k):
        p = []
        while seen[k][1] is not None: p.append(seen[k][2]); k = seen[k][1]
        return p[::-1]
    par = [(seen[k][0], path(k)) for k in seen if seen[k][3] == DEPTH]
    pickle.dump(par, open(A_[3], 'wb')); print('depth-%d parents' % DEPTH, len(par))
elif MODE == 'scan':
    par = pickle.load(open(A_[2], 'rb')); s, e = int(A_[3]), min(int(A_[4]), len(par)); cnt = collections.Counter()
    for i in range(s, e):
        if left() < 15: print('RANGE', s, i, dict(cnt), 'BUDGET spent; resume at', i); sys.exit(2)
        x, p = par[i]; cnt['nodes'] += 1
        for kind, v, st in jcands(x):
            cnt['J:' + st] += 1
            if st == 'kept': print('KEPT index', i, 'path', p, kind, v, 'parent rels', rels(x), flush=True)
    print('RANGE', s, e, dict(cnt), 'DONE %.0fs' % (time.time() - T00))
