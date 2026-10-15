"""Round 066 (skeptic), T5: which side is the repo's End(T)?  H = dim Hom_K(T_i,T_j) (covariant P_x, Hom(P_x,P_y)=paths x~>y,
T_k = P_k -> sum_{k->h} P_h, round-050 complex code); C = ap.cartanMatrix (C[j][i] = #paths i~>j mod I, i.e. dim e_j A e_i).
Part 1: the Oppermann/Ladkani APR example A3/rad^2 (E-181), printed in full, parent and child.
Part 2: sweep all LNAs and duals of length N, every vertex with an out-arrow (no gate): tally child C == H^T / == H / both / neither,
  on steps where H != H^T (so the test has power); and a quiver-level check: every child arrow i->j has H[i][j]>0 (the transposed
  reading would need H[j][i]>0); count arrows violating each reading.  Also the parent: C(parent) == H(A)^T by definition.
Usage: .venv/bin/python workshop/rounds/066/skeptic_side.py N [N ...]"""
import sys
sys.path.insert(0, '.')
from collections import Counter
_st = open('workshop/rounds/050/skeptic_tilt.py').read()
exec(_st[_st.index('from fractions'):_st.index('def stepTest')])
exec(_st[_st.index('def stepTest'):_st.index('def stepTest') + 400].split('\n\nclass')[0])
from quivermutation import nakayama as nk, pathAlgebra

def mat(M): return '\n'.join('   ' + ' '.join('%2d' % x for x in r) for r in M)

def step(alg, v):
    q = alg.quiver; rels = procedure.relationsFrom(alg); V = sorted(q.nodes)
    nq, nr = procedure.mutateAtVertex(q, rels, v)
    C = ap.cartanMatrix(nq, nr, exact=True)
    _, res = stepTest(alg, v, -1)
    return V, nq, nr, C, res

# Part 1
alg = nk.LinearNakayamaAlgebra(3, '2')
V = sorted(alg.quiver.nodes); q = alg.quiver; rels = procedure.relationsFrom(alg)
print('PART 1  A3/rad^2: 1->2->3, relation psi.phi = 0 (path 1,2,3)')
print(' parent arrows', sorted(q.edges(keys=True)), 'rels', rels)
Cp = ap.cartanMatrix(q, rels, exact=True)
print(' Cartan(parent) C[j][i]=#paths i~>j:'); print(mat(Cp))
V, nq, nr, C, res = step(alg, 1)
H = res[0]
print(' H = dim Hom_K(T_i,T_j) for T = (P_1->P_2) + P_2 + P_3 (rows i, cols j):'); print(mat(H))
print(' Hom(T,T[-1]), Hom(T,T[1]) total:', sum(map(sum, res[-1])), sum(map(sum, res[1])))
print(' child arrows', sorted(nq.edges(keys=True)), 'rels', nr)
print(' Cartan(child) C:'); print(mat(C))
n = len(V)
print(' C == H^T:', all(C[i][j] == H[j][i] for i in range(n) for j in range(n)), '  C == H:', all(C[i][j] == H[i][j] for i in range(n) for j in range(n)))
# Part 2
for N in map(int, sys.argv[1:]):
    cnt = Counter(); arrBad = Counter()
    for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
        for tag, a in (('lna', lna), ('dual', pathAlgebra.dualPathAlgebra(lna))):
            Vv = sorted(a.quiver.nodes); n = len(Vv); A0 = Alg(a)
            Tp = complexT(A0, Vv[0], Vv, 0)
            for v in Vv:
                if not ap.arrowsOutOf(a.quiver, v): continue
                try: V, nq, nr, C, res = step(a, v)
                except Exception as e: cnt['error'] += 1; continue
                H = res[0]
                eqT = all(C[i][j] == H[j][i] for i in range(n) for j in range(n))
                eq = all(C[i][j] == H[i][j] for i in range(n) for j in range(n))
                sym = all(H[i][j] == H[j][i] for i in range(n) for j in range(n))
                cnt[('H symmetric' if sym else 'H NOT symmetric', 'C=H^T' if eqT else '', 'C=H' if eq else '', 'neither' if not (eq or eqT) else '')] += 1
                idx = {x: k for k, x in enumerate(V)}
                rq, rr = procedure.reduce(nq, nr)   # arrows removed by a length-1 relation are not arrows of End(T)
                for (t, h, _) in rq.edges(keys=True):
                    arrBad['arrows'] += 1
                    if not H[idx[t]][idx[h]] > 0: arrBad['arrow i->j with H[i][j]=0 (repo reading violated)'] += 1
                    if not H[idx[h]][idx[t]] > 0: arrBad['arrow i->j with H[j][i]=0 (transposed reading violated)'] += 1
            # parent: C(A) == H(A)^T by construction; check at all vertices once per algebra
    print('PART 2 N =', N)
    for k, c in sorted(cnt.items(), key=str): print('  ', c, k)
    for k, c in sorted(arrBad.items()): print('  ', c, k)
