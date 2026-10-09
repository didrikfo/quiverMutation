"""Round 041 (theorist): is the child B of the n = 4 example in the tilting-only class of A?  Bounded BFS from A using only steps that are gate-admitted, legal and have J = 0 (J = 0 <=> Cartan congruence <=> tilting, E-121/E-128).
Prints whether canonicalKey(B) is reached, number of algebras. Usage: theorist_bfs4.py [max]"""
import sys
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); mx = int(ARGV[1]) if len(ARGV) > 1 else 3000
A = build([(1,2),(1,3),(2,3),(2,4),(3,4)], [[[1,2,3,4],[1,3,4]]], [1,2,3,4])
B = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(A, 3))
kB = fingerprint.canonicalKey(B); seen = set(); fr = [A]; n = 0; hit = None; full = True
while fr and n < mx:
    cur, fr[:] = fr[:], []
    for alg in cur:
        k = fingerprint.canonicalKey(alg)
        if k is not None and k in seen: continue
        seen.add(k); n += 1
        if k == kB: hit = n
        if list(nx.simple_cycles(alg.quiver)): continue
        for w in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, w) or perI(alg, w): continue
            raw = mutation.quiverMutationAtVertex(alg, w)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            fr.append(reduction.reducePathAlgebra(raw))
    if n >= mx: full = False
print('tilting-only (J = 0) BFS from A: algebras', n, 'frontier left', len(fr), 'B reached at', hit, 'closed' if not fr else 'NOT closed')
