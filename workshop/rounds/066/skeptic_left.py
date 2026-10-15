"""Round 066 (skeptic), response to referee item 2: check procedure.mutateLeftAtVertex against the dual construction.
Dual construction = the complex T of skeptic_tilt.py built on the OPPOSITE algebra B = A^op (dualPathAlgebra), i.e. Hom(-,A)-side objects of A.
H_B[i][j] = dim Hom_K(T_i,T_j) there.  Expected if mutateLeftAtVertex(A,v) is the dual (right) mutation: Cartan(child_left) == H_B (not its transpose),
and every arrow i->j of the reduced left child has H_B[j][i] > 0 (left child = opposite of B's child; both readings are counted).  Part 1: A3/rad^2 at each vertex with an in-arrow, printed.  Part 2: sweep LNAs n = args.
Same H code as skeptic_side.py (exec of workshop/rounds/050/skeptic_tilt.py, as theorist_enddim.py does).
Usage: .venv/bin/python workshop/rounds/066/skeptic_left.py N [N ...]"""
import sys
sys.path.insert(0, '.')
from collections import Counter
_st = open('workshop/rounds/050/skeptic_tilt.py').read()
exec(_st[_st.index('from fractions'):_st.index('def stepTest')])
exec(_st[_st.index('def stepTest'):_st.index('def stepTest') + 400].split('\n\nclass')[0])
from quivermutation import nakayama as nk, pathAlgebra

def mat(M): return '\n'.join('   ' + ' '.join('%2d' % x for x in r) for r in M)

def left(alg, v):
    q = alg.quiver; rels = procedure.relationsFrom(alg)
    nq, nr = procedure.mutateLeftAtVertex(q, rels, v)
    C = ap.cartanMatrix(nq, nr, exact=True)
    B = pathAlgebra.dualPathAlgebra(alg)
    V, res = stepTest(B, v, -1)
    return V, nq, nr, C, res

def hasIn(q, v): return any(h == v for _, h, _ in q.edges(keys=True))

alg = nk.LinearNakayamaAlgebra(3, '2')
print('PART 1  A3/rad^2 (1->2->3, psi.phi = 0), mutateLeftAtVertex')
for v in (2, 3):
    V, nq, nr, C, res = left(alg, v); H = res[0]; n = len(V)
    print(' v =', v, ' child arrows', sorted(nq.edges(keys=True)), 'rels', nr)
    print('  H_B:'); print(mat(H)); print('  Cartan(left child):'); print(mat(C))
    print('  Hom(T,T[+-1]) total', sum(map(sum, res[-1])), sum(map(sum, res[1])),
          ' C == H_B:', all(C[i][j] == H[i][j] for i in range(n) for j in range(n)),
          ' C == H_B^T:', all(C[i][j] == H[j][i] for i in range(n) for j in range(n)))
for N in map(int, sys.argv[1:]):
    cnt = Counter()
    for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
        for tag, a in (('lna', lna), ('dual', pathAlgebra.dualPathAlgebra(lna))):
            for v in sorted(a.quiver.nodes):
                if not hasIn(a.quiver, v): continue
                try: V, nq, nr, C, res = left(a, v)
                except Exception as e: cnt['error'] += 1; continue
                H = res[0]; n = len(V); idx = {x: k for k, x in enumerate(V)}
                sym = all(H[i][j] == H[j][i] for i in range(n) for j in range(n))
                cnt['steps'] += 1
                cnt['H symmetric (no power)'] += sym
                cnt['C == H_B'] += all(C[i][j] == H[i][j] for i in range(n) for j in range(n))
                cnt['C == H_B^T'] += all(C[i][j] == H[j][i] for i in range(n) for j in range(n))
                cnt['Hom(T,T[+-1]) != 0'] += (sum(map(sum, res[-1])) + sum(map(sum, res[1])) > 0)
                rq, rr = procedure.reduce(nq, nr)
                for (t, h, _) in rq.edges(keys=True):
                    cnt['arrows'] += 1
                    cnt['arrow i->j with H_B[i][j]=0 (same-index reading)'] += not H[idx[t]][idx[h]] > 0
                    cnt['arrow i->j with H_B[j][i]=0 (arrows reversed back, expected reading)'] += not H[idx[h]][idx[t]] > 0
    print('PART 2 N =', N)
    for k, c in cnt.items(): print('  ', c, k)
