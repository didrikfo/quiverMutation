"""Round 137, theorist, T11: the injectivity test (`mutation.isTiltingAtVertex`) equals the Hom criterion on acyclic
quivers, at EVERY vertex (gate-admitted or not).  The proof is in workshop/rounds/137/theorist.md; this checks its three
lemmas separately, per (algebra, vertex, side), with the round-129 Hom engine (chain maps mod homotopy, exact rationals):
    L1  H[1] = H[2] = H[-2] = 0                    (the arrows out of k are an add(D)-approximation; degrees)
    L2  H[-1] = sum_i J_i                          (acyclic: Hom(Y, Y[-1]) = 0)
    THM (H[m] = 0 for m = +-1, +-2)  ==  isTiltingAtVertex
    COR isTiltingAtVertex and k has an out-arrow  =>  gate (mutationIsPossibleAtVertex)
Algebras: round-131 generator (n = 4-7, parallel arrows, zero relations and binomials p + c q, c in {-1, 1, 2}), F side on A,
R side as F on the library's dual.  Power: E-080 (inj NO, Hom NO, gate YES), abd = acd (all YES); planted Y with one arrow
dropped (L1 must fail: H[1] > 0); planted "path-only" injectivity (= the gate's test) must disagree with Hom somewhere.
usage (repo root):  .venv/bin/python -W ignore workshop/rounds/137/theorist_t11.py [controls | sweep N SEED]
"""
import sys, random, collections, warnings
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/129'); sys.path.insert(0, 'workshop/rounds/131')
warnings.simplefilter('ignore')
import theorist_jsilt as J
import theorist_thin as TH
from quivermutation import mutation, pathAlgebra, arrowPaths as ap

def H(A, k, drop=0):
    C = J.Cat(A); T = J.T_plus(A, k, drop)
    return {m: sum(J.homdim(C, X[1], Y[1], m) for X in T for Y in T) for m in (-2, -1, 1, 2)}

def pathOnly(A, k):
    """planted weaker test: no single nonzero PATH into k is killed by every out-arrow (the gate's clause)."""
    C = J.Cat(A); outs = ap.arrowsOutOf(A.quiver, k)
    for i in A.quiver.nodes:
        if i == k: continue
        for p in C.basis(i, k):
            if not any(C.mul({p: 1}, {(b,): 1}, i, k, b[1]) for b in outs): return False
    return True

def one(A, k, tab, bad, tag):
    h = H(A, k); jl = J.J_lib(A, k); inj = mutation.isTiltingAtVertex(A, k)
    gate = mutation.mutationIsPossibleAtVertex(A, k); out = bool(ap.arrowsOutOf(A.quiver, k))
    hom = all(v == 0 for v in h.values())
    L1 = h[1] == h[2] == h[-2] == 0; L2 = h[-1] == jl
    thm = hom == inj; cor = (not (inj and out)) or gate
    tab[(out, gate, inj, hom)] += 1
    tab['pathOnly != Hom'] += (pathOnly(A, k) != hom)
    for name, ok in (('L1', L1), ('L2', L2), ('THM', thm), ('COR', cor)):
        if not ok: bad.append((tag, k, name, h, jl, inj, gate))

def controls():
    sq = [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]
    E080 = J.build(sq, [{(0, 2, 4): 1, (1, 3, 4): -1}]); SQ = J.build(sq, [{(0, 2): 1, (1, 3): -1}])
    for name, A in (('E-080 abde=acde @4', E080), ('abd=acd @4', SQ)):
        print(name, 'H', H(A, 4), 'J', J.J_lib(A, 4), 'inj', mutation.isTiltingAtVertex(A, 4),
              'gate', mutation.mutationIsPossibleAtVertex(A, 4), 'pathOnly', pathOnly(A, 4))
    print('PLANTED drop arrow, E-080 @1: H', H(E080, 1, drop=1), '(L1 must fail)')
    print('PLANTED drop arrow, 1->2->{3,4} @2: H', H(J.build([(1, 2), (2, 3), (2, 4)], []), 2, drop=1), '(L1 must fail)')

def sweep(N, seed):
    rng = random.Random(seed); tab = collections.Counter(); bad = []; cases = 0
    for a in range(N):
        alg = TH.randomAlg(rng, monomial=rng.random() < 0.3)
        P = TH.toLib(alg)
        for side, A in (('F', P), ('R', pathAlgebra.dualPathAlgebra(P))):
            for k in sorted(A.quiver.nodes):
                one(A, k, tab, bad, (seed, a, side)); cases += 1
    print('cases', cases, '(algebras %d, seed %d)' % (N, seed))
    print('(has out-arrow, gate, inj, Hom tilting) -> count')
    for k, v in sorted(((k, v) for k, v in tab.items() if isinstance(k, tuple)), key=str): print('  ', k, v)
    print('planted pathOnly test disagrees with Hom on', tab['pathOnly != Hom'], 'cases')
    print('lemma failures:', len(bad), bad[:10])

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'controls'
    if mode == 'controls': controls()
    else: sweep(int(sys.argv[2]), int(sys.argv[3]))
