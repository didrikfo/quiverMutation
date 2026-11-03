"""Round 111 toolsmith, T10: YES-only search for P in GL_7(Z) with P C_A P^T = C_B for the 89 key-kept J != 0, D = 0 Cartan steps of n = 7 key class 1
(checkpoint /tmp/skeptic_c1_13_ck1.pkl of workshop/rounds/069/skeptic_c1_13.py, 45 000 expansions).  R = P^T lies in the commutant of Phi (R Phi_B = Phi_A R);
if Phi has 7 distinct eigenvalues the commutant has rank 7; R = sum c_i B_i, c in a box; test R^T C_A R == ... exactly.  A hit proves Z-congruence; a miss in the box proves nothing.
Controls: (YES) child of a planted unimodular P0: C' = P0 C_A P0^T must be found; (NO-where-true) the 270 key-changed steps are skipped (Phi differs, no commutant relation).
usage (repo root): .venv/bin/python workshop/rounds/111/toolsmith_c1_13_yes.py [box]"""
import sys, pickle, itertools, collections, random
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
import numpy as np, sympy as sp
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
BOX = int(ARGV[1]) if len(ARGV) > 1 else 2
def phi(C): M = sp.Matrix(C.tolist()); return -(M.T).inv() * M
def commutant(PA, PB):
    """integer basis of {R : R PB = PA R}"""
    n = 7; eqs = []
    unk = sp.symbols('r0:49'); R = sp.Matrix(n, n, unk)
    E = R * PB - PA * R
    A = sp.Matrix([[sp.diff(e, u) for u in unk] for e in E])
    ns = A.nullspace(); out = []
    for v in ns:
        d = sp.ilcm(*[x.q for x in v]); w = [int(x * d) for x in v]; g = sp.igcd(*w); w = [x // g for x in w]
        out.append(np.array(w, dtype=np.int64).reshape(n, n))
    return out
def search(CA, CB, box):
    PA, PB = phi(CA), phi(CB)
    if sp.Poly(PA.charpoly().as_expr()) != sp.Poly(PB.charpoly().as_expr()): return None, 'phi differs'
    Bs = commutant(PA, PB); k = len(Bs)
    if k == 0: return None, 'empty commutant'
    SA = (CA + CA.T).astype(np.int64); SB = (CB + CB.T).astype(np.int64); CAi = CA.astype(np.int64); CBi = CB.astype(np.int64)
    rng = range(-box, box + 1)
    for c in itertools.product(rng, repeat=k):
        if not any(c): continue
        R = sum(ci * B for ci, B in zip(c, Bs))
        P = R.T
        if abs(round(np.linalg.det(P.astype(float)))) != 1: continue
        if (P @ CAi @ P.T == CBi).all(): return P, 'found rank %d' % k
    return None, 'no hit in box %d rank %d' % (box, k)
if __name__ == '__main__':
    seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
    rows = sorted(byCart.values(), key=lambda r: r['nexp']); rows = [r for r in rows if r['kept'] and not r['tp']]
    print('key-kept failing steps', len(rows))
    # control: planted unimodular congruence
    random.seed(1); r0 = rows[0]; CA = r0['CA']
    U = np.eye(7, dtype=np.int64)
    for _ in range(6):
        i, j = random.sample(range(7), 2); U[i] += random.choice([-1, 1]) * U[j]
    Cc = U @ CA @ U.T
    print('control (planted P0 with det', round(np.linalg.det(U)), '):', search(CA, Cc, BOX)[1])
    res = collections.Counter(); misses = []
    for i, r in enumerate(rows):
        P, msg = search(r['CA'], r['CB'], BOX); res[msg.split(' in')[0].split(' rank')[0]] += 1
        if P is None: misses.append((i, r['nexp'], msg))
        if i < 13 or P is None: print(i, r['nexp'], r['v'], msg, flush=True)
    print(dict(res)); print('misses', len(misses))
