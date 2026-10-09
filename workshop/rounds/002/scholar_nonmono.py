"""Round 002 (scholar): negative controls for Ladkani 2.3(c) ("tilt") as coded in
workshop/rounds/001/scholar_h015.py.

Part A  every LNA of length n and its opposite, every vertex with an arrow out:
        gate (mutationIsPossibleAtVertex) x tilt.  Negative control.
Part B  BFS (unguarded, distinct algebras) to --depth; count parents skipped
        (cyclic quiver / baseKey None); on every NON-MONOMIAL parent (some relation
        has >= 2 terms) tabulate gate x tilt x cong over every vertex with an arrow out.
Usage: N=5 DEPTH=3 SHOW=6 .venv/bin/python workshop/rounds/002/scholar_nonmono.py
"""
import argparse, sys, time
from collections import Counter
sys.path.insert(0, '.')
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
sys.argv = ['x']
exec(compile(src, 'h015', 'exec'))   # tiltingPlus, cartan, rplus, imports

def isNonMonomial(alg):
    return any(len(r) >= 2 for r in procedure.relationsFrom(alg))

def run(n, depth, show):
    starts = []
    for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
        starts.append(lna); starts.append(pathAlgebra.dualPathAlgebra(lna))
    A = Counter()
    for s in starts:
        for v in sorted(s.vertices()):
            if not ap.arrowsOutOf(s.quiver, v): continue
            A[(bool(mutation.mutationIsPossibleAtVertex(s, v)), bool(tiltingPlus(s.quiver, procedure.relationsFrom(s), v)))] += 1
    print('PART A n=%d (gate,tilt):' % n, dict(A))
    seen = {}; frontier = []
    def add(alg):
        key = fingerprint.canonicalKey(alg)
        if key is not None:
            if key in seen: return
            seen[key] = 1
        frontier.append(alg)
    for s in starts: add(s)
    T = Counter(); skipped = Counter(); nonmono = 0; expanded = 0; ex = []; t0 = time.time()
    for d in range(depth):
        cur, frontier[:] = frontier[:], []
        for alg in cur:
            expanded += 1
            if search._coxeterKeyOrNone(alg) is None: skipped['baseKey None'] += 1; continue
            if list(nx.simple_cycles(alg.quiver)): skipped['cyclic quiver'] += 1; continue
            nm = isNonMonomial(alg); nonmono += nm
            rels = procedure.relationsFrom(alg); verts = sorted(alg.vertices())
            for v in verts:
                if not ap.arrowsOutOf(alg.quiver, v): continue
                gate = bool(mutation.mutationIsPossibleAtVertex(alg, v))
                tilt = bool(tiltingPlus(alg.quiver, rels, v))
                cong = None
                child = None
                if gate:
                    child = mutation.quiverMutationAtVertex(alg, v)
                    if any(ap.isIllegalRelation(child.quiver, r) for r in procedure.relationsFrom(child)):
                        if nm: T[('nonmono', 'gate', 'illegal-child', tilt)] += 1
                        continue
                    child = reduction.reducePathAlgebra(child)
                    R = rplus(alg, v, verts)
                    cong = bool((R.dot(cartan(alg)).dot(R.T) == cartan(child)).all())
                    add(child)
                if nm:
                    T[('gate' if gate else 'refused', 'tilt' if tilt else 'NOTtilt', cong)] += 1
                    if not tilt and len(ex) < 200: ex.append((d, gate, v, alg.rels, [dict(r) for r in rels if len(r) >= 2][:3]))
        print('depth', d + 1, 'expanded', expanded, 'nonmonomial parents', nonmono, '%.0fs' % (time.time() - t0), flush=True)
    print('skipped parents:', dict(skipped))
    for k, v in sorted(T.items(), key=str): print(k, v)
    print('non-monomial parents with a NOTtilt vertex:', len(ex), ' of which gate-admitted:', sum(1 for e in ex if e[1]))
    for e in ex[:show]: print(e)
    for e in [e for e in ex if e[1]][:show]: print('GATE-ADMITTED NOTtilt', e)

if __name__ == '__main__':
    import os
    run(int(os.environ.get('N', 5)), int(os.environ.get('DEPTH', 3)), int(os.environ.get('SHOW', 6)))
