"""Round 041 (experimentalist): for every gate-admitted step with J_i != 0 (perI(alg,v) nonempty) record whether the CHILD
(reducePathAlgebra of the mutation, key computed on the child by search._coxeterKeyOrNone) has key == parent's key / == class key / is an LNA key at n.
Usage: experimentalist_keyoff.py n cls budget_sec maxdepth [off]
  default (guard ON for expansion): BFS expands only children whose key == class key (as rounds/039/theorist_d3walk.py); every J!=0 step is TESTED, none is filtered.
  'off': expansion also follows children failing the key (up to maxdepth levels, canonicalKey dedup) -> positive-control hunt."""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, cls, budget, maxdepth = int(_a[1]), int(_a[2]), float(_a[3]), int(_a[4]); off = len(_a) > 5 and _a[5] == 'off'
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]; lnakeys = set(classes)
t0 = time.time(); seen = set(); frontier = []; nalg = 0; stop = False; level = 0
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
tab = Counter(); ex = []; nj = 0; dropped = 0; illegal = 0; ngate = 0
while frontier and not stop and level <= maxdepth:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget: stop = True; break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        V = sorted(alg.quiver.nodes); pk = search._coxeterKeyOrNone(alg)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            ngate += 1
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): illegal += 1; continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != V: dropped += 1; continue
            ck = search._coxeterKeyOrNone(ch)
            J = perI(alg, v)
            if J:
                nj += 1
                tab[(pk == base, ck == pk, ck == base, ck in lnakeys)] += 1
                if ck == pk or ck in lnakeys: ex.append((v, J, ck, pk, sorted(alg.quiver.edges(keys=True)), repr(procedure.relationsFrom(alg))))
            if ck == base or off: add(ch)
    level += 1
print('n', n, 'cls', cls, 'seeds', len(classes[base]), 'off' if off else 'on', 'expanded', nalg, 'levels', level, 'stopped', stop, 'gate-admitted', ngate,
      'illegal', illegal, 'vertexset-dropped', dropped, 'J!=0 steps', nj, '%.0fs' % (time.time()-t0))
print('(parent key==class, child==parent, child==class, child is an LNA key at n): count')
for k, c in sorted(tab.items(), key=str): print('  ', k, c)
for e in ex[:6]: print('EX', e)
