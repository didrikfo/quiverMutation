"""Cartan congruence check for replayed parents of E-084 (round 015, theorist).
Usage: theorist_cartan.py n class     -- replays every recorded line of workshop/rounds/014/scholar_walk_n{n}_c{class}.txt
 ('REJECTIONS' lines: gate True, tiltingPlus False;  'M' lines: gate True, tiltingPlus True, key moved).
For each parent: gate, tiltingPlus, cong = (C(child) == R C(parent) R^T) for the reduced child AND the unreduced child,
key of child == key of parent, number of parallel arrows / duplicated relation terms in the parent."""
import sys, ast
n, c = int(sys.argv[1]), int(sys.argv[2])
lines = []
mode = None
for l in open('workshop/rounds/014/scholar_walk_n%d_c%d.txt' % (n, c)):
    if l.startswith('REJECTIONS'): mode = 'R'
    elif l.startswith('GATE+TILT'): mode = 'M'
    elif mode == 'R' and l.startswith('('): lines.append(('R', ast.literal_eval(l.strip())))
    elif mode == 'M' and l.startswith('M '): lines.append(('M', ast.literal_eval(l[2:].strip())))
sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[c]
for kind, t in lines:
    if kind == 'R': d, rels, v, g, path = t
    else: d, rels, v, path = t
    alg = classes[base][path[0]]
    for w in path[1:]:
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, w))
    pk = search._coxeterKeyOrNone(alg)
    verts = sorted(alg.vertices())
    rl = procedure.relationsFrom(alg)
    gate = mutation.mutationIsPossibleAtVertex(alg, v)
    tp = tiltingPlus(alg.quiver, rl, v)
    raw = mutation.quiverMutationAtVertex(alg, v)
    ch = reduction.reducePathAlgebra(raw)
    R = rplus(alg, v, verts); C = cartan(alg)
    RCR = R.dot(C).dot(R.T)
    congRed = bool((RCR == cartan(ch)).all()) if len(ch.vertices()) == len(verts) else 'dim'
    congRaw = bool((RCR == cartan(raw)).all()) if len(raw.vertices()) == len(verts) else 'dim'
    edges = [(a, b) for a, b, *_ in alg.quiver.edges(keys=True)]
    par = len(edges) - len(set(edges))
    # sign/size of the discrepancy
    diff = (RCR - cartan(ch)) if congRed != 'dim' else None
    print(kind, 'depth', d, 'v', v, 'gate', gate, 'tiltingPlus', tp, 'congReduced', congRed, 'congUnreduced', congRaw,
          'keyEq', search._coxeterKeyOrNone(ch) == pk, 'parallel', par, 'outdeg(v)', len(ap.arrowsOutOf(alg.quiver, v)),
          'entries<0 in C(child)', int((cartan(ch) < 0).sum()), 'diff', None if diff is None else sorted(set(int(x) for x in diff.flatten() if x)))
    print('   parent rels', alg.rels)
