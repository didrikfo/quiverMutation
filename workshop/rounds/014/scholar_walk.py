"""T5 round 014: guarded BFS from the LNAs (and their relation duals) of length n, split by Coxeter key class.
At every distinct algebra reached and every gate-admitted vertex (mutationIsPossibleAtVertex) record whether
tiltingPlus (= AI 2.32(b) = Ladkani 2.3(c), from workshop/rounds/001/scholar_h015.py) holds; a gate-admitted
step with tiltingPlus False is a REJECTION (the thing asked for). The walk follows only guard-passing steps
(child Coxeter key = parent key), so it is a guarded walk; each key class is independent.
  scholar_walk.py n --plan                       class sizes (starts per key), no walking
  scholar_walk.py n --class I [--depth D] [--budget-sec S]   walk one class (I = index in --plan order)
Exit 0 if closed (frontier empty), 2 if stopped by depth or budget (counts then partial)."""
import argparse, sys, time
from collections import Counter
sys.path.insert(0, '.'); sys.argv = sys.argv[:1] + ['0'] if False else sys.argv
_argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
ap_ = argparse.ArgumentParser()
ap_.add_argument('n', type=int); ap_.add_argument('--plan', action='store_true')
ap_.add_argument('--class', type=int, dest='cls', default=-1); ap_.add_argument('--depth', type=int, default=99)
ap_.add_argument('--stop-on-reject', action='store_true', dest='sor')
ap_.add_argument('--budget-sec', type=float, default=0, dest='budget'); ap_.add_argument('--show', type=int, default=5)
a = ap_.parse_args()

classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(a.n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        k = search._coxeterKeyOrNone(alg)
        classes.setdefault(k, []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
if a.plan:
    for i, k in enumerate(order): print(i, len(classes[k]), k)
    sys.exit(0)
starts = classes[order[a.cls]]; base = order[a.cls]
seen = set(); frontier = []; paths = {}
def add(alg, path=()):
    key = fingerprint.canonicalKey(alg)
    if key is not None:
        if key in seen: return False
        seen.add(key)
    frontier.append(alg); paths[id(alg)] = path; return True
for i, s in enumerate(starts): add(s, (i,))
tab = Counter(); rej = []; mism = []; t0 = time.time(); expanded = 0; stop = False; nonmono = 0
for d in range(a.depth):
    cur, frontier[:] = frontier[:], []
    if not cur: break
    for alg in cur:
        if a.budget and time.time() - t0 > a.budget: stop = True; break
        expanded += 1
        if list(nx.simple_cycles(alg.quiver)): tab['parent-cyclic'] += 1; continue
        rels = procedure.relationsFrom(alg)
        if any(len(r) > 1 for r in rels): nonmono += 1
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            t = tiltingPlus(alg.quiver, rels, v)
            child = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(child.quiver, r) for r in procedure.relationsFrom(child)):
                tab[('illegal', t)] += 1; continue
            child = reduction.reducePathAlgebra(child)
            ck = search._coxeterKeyOrNone(child)
            guard = ck == base
            tab[('guard' if guard else 'noguard', 'tilt' if t else 'NOTtilt')] += 1
            if t and not guard: mism.append((d, alg.rels, v, paths[id(alg)]))
            if not t:
                rej.append((d, alg.rels, v, guard, paths[id(alg)]))
            if guard: add(child, paths[id(alg)] + (v,))
    print('depth', d + 1, 'expanded', expanded, 'next', len(frontier), 'nonmono parents', nonmono, '%.0fs' % (time.time() - t0), flush=True)
    if stop or (a.sor and rej): break
stop = stop or bool(a.sor and rej)
closed = not frontier and not stop
print('n', a.n, 'class', a.cls, 'key', base, 'starts', len(starts), 'distinct algebras', len(seen), 'closed', closed)
for k, v in sorted(tab.items(), key=str): print(k, v)
print('REJECTIONS (gate-admitted, tiltingPlus False):', len(rej))
for r in rej[:a.show]: print(r)
print('GATE+TILT BUT KEY MOVES:', len(mism))
for r in mism[:a.show]: print('M', r)
sys.exit(0 if closed else 2)
