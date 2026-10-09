"""Self-test of skeptic_tilt.py: on A itself (T = A) H must equal the Cartan matrix; on LNA steps known to be tilting."""
import sys; sys.path.insert(0, '.')
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
from quivermutation import nakayama as nk, invariants, pathAlgebra
import numpy as np
S = -1
for lna in list(nk.LinearNakayamaAlgebra.allOfLength(5))[:6]:
    A = Alg(lna); V = sorted(lna.quiver.nodes)
    T = {i: ({0: [i]}, {}) for i in V}
    H = [[homdim(A, T[i], T[j], 0) for j in V] for i in V]
    C = np.array(invariants.cartanMatrix(lna, exact=True).tolist(), dtype=int)
    print('Hom(A,A)==Cartan', (np.array(H) == C).all(), (np.array(H) == C.T).all())
    for k in V:
        V_, r = stepTest(lna, k, S)
        print(' k', k, 'H-1 zero', not np.array(r[-1]).any(), 'H+1 zero', not np.array(r[1]).any())
