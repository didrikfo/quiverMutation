"""T10(i): one-step reverse test. For each failing E-149 child B (from a skeptic_collect pickle), mutate B at the same vertex v
(plain, and via the opposite algebra) and report gate / J / key / whether the result is the parent (canonicalKey) or parent's opposite.
If a J = 0 gate-admitted reverse step returned the parent, B -> P would be a tilting step and B would be in the parent's class.
usage: theorist_reverse.py in.pkl   (run from repo root)"""
import sys, pickle, ast
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
def mk(edges, reprr, nodes):
    arrows = [(a, b) for a, b, k in edges]; rels = []
    for d in ast.literal_eval(reprr): rels.append([[p[0][0]] + [a[1] for a in p] for p in d])
    return build(arrows, rels, nodes)
recs = [x for x in pickle.load(open(ARGV[1], 'rb')) if x['kind'] == 'fail']
for i, r in enumerate(recs):
    A = mk(r['pe'], r['pr'], r['V']); v = r['v']
    B = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(A, v))
    kA = fingerprint.canonicalKey(A); kAo = fingerprint.canonicalKey(pathAlgebra.dualPathAlgebra(A))
    out = []
    for name, X in (('same', B), ('opp', pathAlgebra.dualPathAlgebra(B))):
        try:
            ok = mutation.mutationIsPossibleAtVertex(X, v)
            if not ok: out.append((name, 'gate-closed')); continue
            Y = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(X, v))
            kY = fingerprint.canonicalKey(Y)
            out.append((name, 'gate-open', 'J', perI(X, v), 'J_is_zero', not perI(X, v), 'is_parent', kY == kA, 'is_opp_parent', kY == kAo))
        except Exception as e:
            out.append((name, 'error', repr(e)[:40]))
    print(i, 'depth', r['depth'], 'v', v, out, flush=True)
