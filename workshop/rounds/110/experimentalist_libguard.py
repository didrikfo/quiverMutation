"""Round 110 experimentalist, T10 / H-015: key-refusal census.  Guarded BFS (expansion only into key-kept children, as the library guard does)
from the LNAs and duals of key class CLS at length N, canonicalKey dedup, wall budget.  EVERY gate-admitted legal vertex-preserving step is classified:
  J = perI(alg, v) nonempty (J != 0) or empty (J = 0);  key moved (the guard REFUSES it) or kept;  tiltingPlus on the parent (Ladkani 2.3(c)).
Prints a table (J, moved, tilting) -> count, and the first refused steps with parent depth, vertex, keys.  A refused step is certainly not a derived
equivalence (key = Coxeter polynomial, invariant of derived equivalence for finite global dimension), so a refusal is always a correct NO.
Usage (repo root): experimentalist_refusal.py N CLS budget_sec [nshow]"""
import sys, time, copy
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
n, cls, budget = int(_a[1]), int(_a[2]), float(_a[3]); nshow = int(_a[4]) if len(_a) > 4 else 3; NREF = 25; cap_ref = []; cap_kept = []
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
t0 = time.time(); seen = set(); frontier = []; tab = Counter(); ex = []; nexp = 0; level = 0; stop = False
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
while frontier and not stop:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget or len(cap_ref) >= NREF: stop = True; break
        nexp += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        V = sorted(alg.quiver.nodes); pk = search._coxeterKeyOrNone(alg)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): tab['illegal'] += 1; continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != V: tab['vertexset-dropped'] += 1; continue
            ck = search._coxeterKeyOrNone(ch)
            if ck is None: tab['cyclic-child'] += 1; continue
            J = bool(perI(alg, v)); moved = ck != pk
            tp = bool(tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v))
            tab[('J!=0' if J else 'J=0', 'REFUSED' if moved else 'kept', 'tilting' if tp else 'NOT-tilting')] += 1
            if moved and len(ex) < nshow: ex.append((level, v, pk, ck, sorted(alg.quiver.edges(keys=True)), repr(procedure.relationsFrom(alg))))
            if moved and len(cap_ref) < NREF: cap_ref.append((copy.deepcopy(alg), v))
            if not moved and tab[('J=0','kept','tilting')] % 700 == 0 and len(cap_kept) < 60: cap_kept.append((copy.deepcopy(alg), v))
            if not moved: add(ch)
    level += 1
print('n', n, 'cls', cls, 'classsize', len(classes[base]), 'classkey', base, 'expanded', nexp, 'levels', level, 'budget-hit', stop, '%.0fs' % (time.time() - t0))
for k, c in sorted(tab.items(), key=str): print('  ', k, c)
for e in ex: print('EX level, v, parentkey, childkey, arrows, rels:', e)

def libvisits(alg, guard):
    seen = []
    search.mutationSearchDepthFirst(copy.deepcopy(alg), 1, printOutput=False, coxeterGuard=guard, visitor=lambda a, mv: seen.append(tuple(mv)))
    return seen
for name, cap, expect in (('refused', cap_ref, False), ('kept', cap_kept, True)):
    agree = 0; offok = 0
    for alg, v in cap:
        g = (v,) in set(libvisits(alg, True)); u = (v,) in set(libvisits(alg, False))
        agree += (g == expect); offok += u
    print('LIBRARY GUARD check', name, 'steps', len(cap), 'library-guarded search visits child =', expect, 'in', agree, '; unguarded search visits child in', offok)
