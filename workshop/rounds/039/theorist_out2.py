"""Round 039 (theorist): hand algebra with J_i != 0, d_i = 3, out(i) = 2, gate-admitted (so out(i) = 3 is not forced by the gate);
test its Coxeter key vs the n=7 LNA keys, the child key, and E-136's C_B = r C_A r^T + H.  Usage: theorist_out2.py"""
import sys
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
_a = ARGV; exec(compile(src, 'bd', 'exec'))
# vertices 1=i 2=a 3=b 4=c 5=d 6=v 7=w
arr = [(1,2),(1,3),(2,4),(2,5),(4,6),(5,6),(3,6),(6,7)]
p1 = [1,2,4,6,7]; p2 = [1,2,5,6,7]; p3 = [1,3,6,7]
cases = {'out2 comm p1=p3 (J=p1-p3)': [[p1, p3]], 'out2 comm p1=p3 and p2=p3 (J dim2)': [[p1,p3],[p2,p3]]}
ks = lnaKeys(7)
for name, rl in cases.items():
    A = build(arr, rl, list(range(1, 8))); v = 6
    print(name, '| gate', mutation.mutationIsPossibleAtVertex(A, v), '| (d,J) at i=1', dims(A, v) if 'dims' in globals() else '')
    J = perI(A, v); Q = A.quiver
    P = ap.allPathsBetween(Q, 1, 6); rels = procedure.relationsFrom(A)
    print('  d_1 =', len(P) - len(ap.idealBasis(Q, rels, 1, 6)), 'J =', J, 'out(1) =', len(ap.arrowsOutOf(Q, 1)))
    print('  key', search._coxeterKeyOrNone(A), 'in LNA keys:', search._coxeterKeyOrNone(A) in ks)
    if mutation.mutationIsPossibleAtVertex(A, v):
        raw = mutation.quiverMutationAtVertex(A, v)
        if not any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)):
            ch = reduction.reducePathAlgebra(raw); kc = search._coxeterKeyOrNone(ch)
            print('  child key == parent key:', kc == search._coxeterKeyOrNone(A))
