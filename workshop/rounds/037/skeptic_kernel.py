"""Round 037 (skeptic): the kernel element J_i of each E-129 row as a combination of path classes e_iAe_v, and which out-arrows it uses.
Usage: skeptic_kernel.py skeptic_rows.pkl"""
import sys, pickle
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/037/skeptic_orbits.py').read().split("algs = []")[0]
sys.argv = ARGV; exec(compile(src, 'orb', 'exec'))
import sympy
for (nalg, v, e, r, J) in pickle.load(open(ARGV[1], 'rb')):
    A = build(e, r); Q = A.quiver; rels = procedure.relationsFrom(A); outs = ap.arrowsOutOf(Q, v)
    print('ROW', nalg, 'v', v, 'outs', outs)
    for i in sorted(J):
        P = ap.allPathsBetween(Q, i, v); piv = ap.idealBasis(Q, rels, i, v)
        # reduce each path mod ideal -> vectors; kernel of (path -> images) modulo ideal: use kernel of stacked [img rows] and report those with nonzero class
        rws = []; keys = {}
        for q in P:
            row = {}
            for b in outs:
                for kk, x in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(Q, rels, i, b[1])).items(): row[(b, kk)] = x
            rws.append(row); [keys.setdefault(k, len(keys)) for k in row]
        M = sympy.Matrix([[sympy.Rational(r.get(k, 0)) for k in keys] for r in rws]) if keys else sympy.zeros(len(P), 0)
        ns = M.T.nullspace()
        print(' i', i, 'paths', [[a[1] for a in q] for q in P], 'ideal rank', len(piv))
        for vec in ns: print('   null (path coefficients):', list(vec), ' image arrows used:', sorted({b for r_ in rws for (b, kk) in r_}))
