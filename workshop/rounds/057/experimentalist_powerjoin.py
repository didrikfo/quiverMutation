"""Round 057, T10 power control for the J = 0 join test (specificity).
Walks from LNA rows A and B (n = 9 default) in two modes, both with the gate and the Coxeter-key guard:
  J0  : every step must also satisfy tiltingPlus (J = 0), the E-155/E-161/E-167 test;
  ALL : gate + key guard only (the search's own walk, J != 0 steps allowed).
Both with right steps from the LNA and from its dual (left steps).  Then a relabelling-aware join (WL hash + VF2, as amerge) of A's reach set with B's.
A certified-equivalent pair (same quipu, F-047) must join (sensitivity); a certified-inequivalent pair with equal Coxeter polynomial
(F-010 pair: Z-conjugacy profile of F-047 differs, Ladkani Cor 3.15) must not join if the J0 test is sound.
Usage (repo root): timeout 10m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py N ROWA ROWB DEPTH [modes=J0,ALL]
   e.g. ... 9 3060000 3304000 4"""
import sys, copy, time, importlib.util
sys.path.insert(0, '.')
spec = importlib.util.spec_from_file_location('am', 'workshop/rounds/054/experimentalist_amerge.py')
am = importlib.util.module_from_spec(spec); spec.loader.exec_module(am)
from quivermutation import nakayama as nk, mutation, procedure, reduction, search, pathAlgebra, lnaMoves, arrowPaths
TP = am.loadTP()

def walk(alg0, depth, j0only, dual=False):
    """{labelled quiverKey: (path, nJ!=0 steps)} reached; path steps are original vertex numbers (negated for left steps)."""
    base = search._coxeterKeyOrNone(alg0)
    found = {}; seen = {}
    start = pathAlgebra.dualPathAlgebra(alg0) if dual else alg0
    def rec(a, path, bad):
        k = search.quiverKey(a if not dual else pathAlgebra.dualPathAlgebra(a))
        if k is not None and (k not in found or len(path) < len(found[k][0])): found[k] = (list(path), bad)
    def go(a, d, path, bad):
        rec(a, path, bad)
        if d == 0 or search.hasOrientedCycle(a): return
        k = search.quiverKey(a)
        if k is not None:
            if seen.get(k, -1) >= d: return
            seen[k] = d
        for v in sorted(a.quiver.nodes):
            if not mutation.mutationIsPossibleAtVertex(a, v): continue
            j0 = bool(TP(a.quiver, procedure.relationsFrom(a), v))
            if j0only and not j0: continue
            m = mutation.quiverMutationAtVertex(copy.deepcopy(a), v)
            if any(arrowPaths.isIllegalRelation(m.quiver, r) for r in procedure.relationsFrom(m)): continue
            m = reduction.reducePathAlgebra(m)
            mk = search._coxeterKeyOrNone(m)
            if mk is not None and mk != base: continue
            go(m, d - 1, path + [(-v if dual else v)], bad + (0 if j0 else 1))
    go(copy.deepcopy(start), depth, [], 0)
    return found

def reachSet(row, n, depth, j0only):
    alg = nk.LinearNakayamaAlgebra(n, row)
    f = walk(alg, depth, j0only, False)
    for k, v in walk(alg, depth, j0only, True).items():
        if k not in f or len(v[0]) < len(f[k][0]): f[k] = v
    return f

if __name__ == '__main__':
    n = int(sys.argv[1]); ra = [int(c) for c in sys.argv[2]]; rb = [int(c) for c in sys.argv[3]]; D = int(sys.argv[4])
    modes = sys.argv[5].split(',') if len(sys.argv) > 5 else ['J0', 'ALL']
    print('pair', sys.argv[2], sys.argv[3], 'n', n, 'depth each side', D)
    for mode in modes:
        t0 = time.time()
        FA = reachSet(ra, n, D, mode == 'J0'); FB = reachSet(rb, n, D, mode == 'J0')
        A = {k: v[0] for k, v in FA.items()}; B = {k: v[0] for k, v in FB.items()}
        hits = am.join(A, {'B': B})
        lnaA = sum(1 for k in A if lnaMoves.asRelLengths(None, n) is None) if False else None
        print(mode, 'reach A', len(A), 'reach B', len(B), 'iso meetings', len(hits),
              'shortest total', hits[0][0] if hits else None, 'sec %.0f' % (time.time() - t0), flush=True)
        for h in hits[:2]: print('   fwd', h[2], 'back', h[3], 'inverse-back', h[4])
        # was any J != 0 step on the winning route?
        if hits: print('   J!=0 steps on the A-side path', FA[hits[0][5]][1], ' B-side', FB[hits[0][6]][1] if hits[0][6] in FB else '?')
