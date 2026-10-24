"""Round 085 (toolsmith; copy of 083 with KEYROOT and DEPTH): , T10: scan the depth-5 parents of the 3060000 J = 0 ball (n = 9) for key-kept perI (J != 0) steps.
usage: toolsmith_d5scan.py build ROOT OUT.pkl         BFS to depth 5 (as skeptic_kept.py find); saves depth-5 nodes (the parents) with parent pointers
       toolsmith_d5scan.py scan IN.pkl START END     jcands-style scan of parents[START:END]: statuses of every gate-admissible perI step; prints KEPT ones
       toolsmith_d5scan.py ctl IN.pkl                positive control: rerun the scan on the one known depth-4 kept parent (E-233) -> must print KEPT R 8
Run from the repository root."""
import sys, time, pickle, collections
A_ = sys.argv; MODE = A_[1]
sys.path.insert(0, '.'); sys.argv = ['x', 'none', '2', '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import os
KEYROOT = os.environ.get('KEYROOT', '3060000'); DEPTH = int(os.environ.get('DEPTH', '5'))
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
if MODE == 'build':
    ROOT, OUT = A_[2], A_[3]; root = nk.LinearNakayamaAlgebra.fromClassName(ROOT)
    seen = {fingerprint.canonicalKey(root): (root, None, None, 0)}; front = [root]; T0 = time.time()
    for d in range(DEPTH):
        nxt = []
        for x in front:
            kx = fingerprint.canonicalKey(x)
            for kind, v, c in moves(x):
                k = fingerprint.canonicalKey(c)
                if k is not None and k not in seen: seen[k] = (c, kx, (kind, v), d + 1); nxt.append(c)
        print('depth', d + 1, 'ball', len(seen), '%.0fs' % (time.time() - T0), flush=True); front = nxt
    def path(k):
        p = []
        while seen[k][1] is not None: p.append(seen[k][2]); k = seen[k][1]
        return p[::-1]
    par = [(seen[k][0], path(k)) for k in seen if seen[k][3] == DEPTH]
    pickle.dump(par, open(OUT, 'wb')); print('depth-%d parents' % DEPTH, len(par))
elif MODE == 'scan':
    par = pickle.load(open(A_[2], 'rb')); s, e = int(A_[3]), int(A_[4]); cnt = collections.Counter(); T0 = time.time()
    for i in range(s, min(e, len(par))):
        x, p = par[i]; cnt['nodes'] += 1
        for kind, v, st in jcands(x):
            cnt['J:' + st] += 1
            if st == 'kept': print('KEPT index', i, 'path', p, kind, v, 'parent rels', rels(x), flush=True)
    print('RANGE', s, e, dict(cnt), '%.0fs' % (time.time() - T0))
elif MODE == 'ctl':
    f = pickle.load(open(A_[2], 'rb'))[0]; print('ctl', jcands(f['parent']))
