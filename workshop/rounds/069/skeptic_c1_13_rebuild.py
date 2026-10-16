"""Round 069 (skeptic): rebuild by hand (arrows + relations only) the first 13 D == 0 class-1 steps recorded by skeptic_c1_13.py
(/tmp/skeptic_c1_13_cls1.pkl, E147 list in discovery order), parent and child; parent without parallel arrows only (the hand constructor
build() cannot take them).  Checks gate, tiltingPlus, key(parent) == key(child), canonicalKey(step child) == canonicalKey(hand child)."""
import sys, pickle, ast
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
def mk(edges, reprr, nodes):
    arrows = [(a, b) for a, b, k in edges]; rels = []
    for d in ast.literal_eval(reprr): rels.append([[p[0][0]] + [a[1] for a in p] for p in d])
    return build(arrows, rels, nodes)
recs = [r for r in pickle.load(open(ARGV[1], 'rb')) if r['tag'] == 'E147'][:13]
for i, r in enumerate(recs):
    if r['par']: print(i, 'parent has parallel arrows', r['par'], '-- not hand-rebuilt'); continue
    N = r['V']; v = r['v']; A = mk(r['pe'], r['pr'], N); B = mk(r['ce'], r['cr'], N)
    ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(A, v))
    print(i, 'v', v, 'arrows', [e[:2] for e in r['pe']], 'gate', mutation.mutationIsPossibleAtVertex(A, v), 'tiltingPlus', bool(tiltingPlus(A.quiver, procedure.relationsFrom(A), v)),
          'key same', search._coxeterKeyOrNone(A) == search._coxeterKeyOrNone(ch) == search._coxeterKeyOrNone(B), 'child==hand child', fingerprint.canonicalKey(ch) == fingerprint.canonicalKey(B), 'J', perI(A, v))
