"""Round 081 (experimentalist), S-1/T10: reconcile E-219 'J != 0 candidates rejected in the balls (0.6%)' with E-175/E-204 'no J != 0 step enters the balls'.
For each node expanded in the key-kept J = 0 ball of ROOT (depth <= DEPTH-1 levels expanded), every gate-admissible (kind, v) with J != 0 (perI) is also MUTATED
and its child's Coxeter key compared with the class key.  usage: experimentalist_jkey.py ROOT DEPTH BUDGET_S   (run from the repository root)
Counts: Jcand (E-219 count), Jkeykept (child key == class key: the T10-relevant 'J != 0 step'), Jkeymoved, Jbad (illegal / vertex set changed / None)."""
import sys, time, collections
A_ = sys.argv; ROOT, DEPTH, BUD = A_[1], int(A_[2]), float(A_[3])
sys.path.insert(0, '.'); sys.argv = ['x', 'none', '2', '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
base = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3060000'))
assert base == search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3304000'))
T0 = time.time()
def jcands(alg):
    """children of gate-admissible J != 0 candidates: list of (kind, v, status)"""
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
                if k == base: print('KEPT J!=0 step: parent', [list(r) for r in procedure.relationsFrom(alg)], 'kind', kind, 'v', v, 'child rels', [list(r) for r in procedure.relationsFrom(c)], 'child key in ball', fingerprint.canonicalKey(c) in seen, 'child tiltingPlus-J ok', flush=True)
            except Exception as e: out.append((kind, v, 'err ' + type(e).__name__))
    return out
root = nk.LinearNakayamaAlgebra.fromClassName(ROOT)
seen = {fingerprint.canonicalKey(root): root}; front = [root]; tot = collections.Counter(); nexp = 0; perlevel = []
for d in range(DEPTH):
    nxt = []; lv = collections.Counter(); done = True
    for x in front:
        if time.time() - T0 > BUD: done = False; break
        nexp += 1; lv['nodes'] += 1
        for kind, v, st in jcands(x): lv['J:' + st] += 1
        for kind, v, c in moves(x):
            k = fingerprint.canonicalKey(c)
            if k is not None and k not in seen: seen[k] = c; nxt.append(c)
    perlevel.append((d, dict(lv), done)); tot.update(lv); print('level', d, dict(lv), 'complete' if done else 'BUDGET-STOP', 'ball', len(seen), '%.0fs' % (time.time() - T0), flush=True)
    if not done: break
    front = nxt
print('TOTAL', dict(tot), 'gate/J/tp/key/ok', stats)
