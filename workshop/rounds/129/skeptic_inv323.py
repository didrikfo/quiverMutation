"""Round 129 (skeptic), T10 blocking check for E-323, part (2)+(3): invariants of the pair 0040033030 / 0400330030 other than SNF(S),
computed without any Smith-normal-form routine, with YES controls whose derived equivalence is certified at Hom level.
Phi := -C^{-1} C^T (C unimodular).  C^{-1} S = 1 - Phi for S = C + C^T, so SNF(S) = coker(Phi - 1); every rank_p f(Phi) below is invariant
under P Phi_A = Phi_B P, P in GL_n(Z) -- the Phi half of H-024's necessity step (Ladkani 1001.4765 Lemma 3.5 + Aihara-Iyama Thm 2.27, E-308),
which needs no lattice theory (Nikulin 1.13.3 is not used).  Invariant vector: for f in FAM and p in PRIMES, rank of f(Phi) over F_p (own
Gaussian elimination).  'x-1' is SNF(S) restated; all other f are the "other" invariants; we report which of them separate the pair.
Walk: the 5 steps of E-323 (F3 R9 F3, F3 (J != 0), and the J = 0 key-kept step that reaches 0400330030, found here), each with the Hom
test of skeptic_hom323.py and End(T)'s Cartan matrix from Hom dims.  Controls: YES = each Hom-tilting step (End(T) vs parent, truth YES by
Rickard); YES = r C r^T for random r in GL_n(Z) (code); planted error = r with det 2 or 3 (caught only by the unimodularity check: no power); E-318's n = 11 pair.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/129/skeptic_inv323.py"""
import sys, random, copy, time
sys.path.insert(0, 'workshop/rounds/129')
import skeptic_hom323 as h
from fractions import Fraction
PRIMES = (2, 3, 5, 7, 11, 13)
FAM = {'x-%d' % a if a >= 0 else 'x+%d' % -a: [-a, 1] for a in range(-4, 5)}          # coefficient lists, low degree first
FAM.update({'x^%d-1' % k: [-1] + [0] * (k - 1) + [1] for k in range(2, 7)})
FAM.update({'x^%d+1' % k: [1] + [0] * (k - 1) + [1] for k in range(2, 7)})
FAM.update({'x^2+x+1': [1, 1, 1], 'x^2-x+1': [1, -1, 1]})

def mul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def inv_int(C):                                   # exact inverse, must be integral (det C = +-1)
    n = len(C); M = [[Fraction(C[i][j]) for j in range(n)] + [Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        M[c] = [v / M[c][c] for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0: M[r] = [a - M[r][c] * b for a, b in zip(M[r], M[c])]
    out = [[M[i][n + j] for j in range(n)] for i in range(n)]
    assert all(v.denominator == 1 for r in out for v in r), 'C not unimodular'
    return [[int(v) for v in r] for r in out]
def phi(C):
    Ct = [list(r) for r in zip(*C)]; X = mul(inv_int(C), Ct); return [[-v for v in r] for r in X]
def poly(Phi, co):
    n = len(Phi); R = [[0] * n for _ in range(n)]; Pk = [[int(i == j) for j in range(n)] for i in range(n)]
    for c in co:
        R = [[R[i][j] + c * Pk[i][j] for j in range(n)] for i in range(n)]; Pk = mul(Pk, Phi)
    return R
def rank_p(M, p):
    M = [[v % p for v in r] for r in M]; rk = 0; n, m = len(M), len(M[0])
    for c in range(m):
        piv = next((r for r in range(rk, n) if M[r][c]), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; iv = pow(M[rk][c], p - 2, p)
        M[rk] = [v * iv % p for v in M[rk]]
        for r in range(n):
            if r != rk and M[r][c]: f = M[r][c]; M[r] = [(a - f * b) % p for a, b in zip(M[r], M[rk])]
        rk += 1
    return rk
def charpoly_mod(Phi, p=1000003):                 # fingerprint of the characteristic polynomial: det(t - Phi) at t = 2..6 mod p
    n = len(Phi); out = []
    for t in range(2, 7):
        M = [[((t if i == j else 0) - Phi[i][j]) % p for j in range(n)] for i in range(n)]; d = 1
        for c in range(n):
            piv = next((r for r in range(c, n) if M[r][c]), None)
            if piv is None: d = 0; break
            if piv != c: M[c], M[piv] = M[piv], M[c]; d = -d
            d = d * M[c][c] % p; iv = pow(M[c][c], p - 2, p)
            for r in range(c + 1, n):
                f = M[r][c] * iv % p; M[r] = [(a - f * b) % p for a, b in zip(M[r], M[c])]
        out.append(d % p)
    return tuple(out)
def invec(C):
    F = phi(C); return {'charpoly': charpoly_mod(F), **{(f, p): rank_p(poly(F, co), p) for f, co in FAM.items() for p in PRIMES}}
def diff(a, b): return sorted(k for k in a if a[k] != b[k])

def lnaC(name):
    N = len(name) + 2; return [list(map(int, r)) for r in h.ct.lnaCartanMatrix(N, tuple(int(c) for c in name))]
def cartan(alg): return [list(map(int, r)) for r in h.invariants.integerCartanMatrix(alg)]
def unimod(n, rng, det=1):
    R = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(4 * n):
        i, j = rng.sample(range(n), 2); c = rng.choice((-1, 1)); R[i] = [a + c * b for a, b in zip(R[i], R[j])]
    if det != 1: R[0] = [det * v for v in R[0]]
    rng.shuffle(R); return R

if __name__ == '__main__':
    t0 = time.time(); N = 12; rng = random.Random(129)
    A, B = '0040033030', '0400330030'
    CA, CB = lnaC(A), lnaC(B)
    vA, vB = invec(CA), invec(CB)
    print('pair %s / %s: brute-force Cartan == table: %s %s' % (A, B, cartan(h.nk.LinearNakayamaAlgebra.fromClassName(A)) == CA,
          cartan(h.nk.LinearNakayamaAlgebra.fromClassName(B)) == CB))
    d = diff(vA, vB)
    print('charpoly(Phi) equal:', vA['charpoly'] == vB['charpoly'])
    print('separating (f, p):', d)
    print('separating, other than f = x-1 (= SNF(S)):', [k for k in d if k[0] != 'x-1'])
    for k in d: print('   %-10s p=%-2d rank %s vs %s' % (k[0], k[1], vA[k], vB[k]))
    # the walk, step by step, with Hom checks and End(T) Cartan from Hom dims
    print('\nE-323 walk from %s (each step: J, Hom-tilting, invariants of End(T) == parent?):' % A)
    x = h.nk.LinearNakayamaAlgebra.fromClassName(A); pk0 = h.search._coxeterKeyOrNone(x); cur = CA
    for mv in ('F3', 'R9', 'F3', 'F3'):
        res = h.check(x, mv, 'walk'); x, _ = h.step(x, mv); C = cartan(x)
        print('      invariants(child) == invariants(parent): %s; == invariants(root): %s; differs at %d of %d (f,p)' %
              (invec(C) == invec(cur), invec(C) == vA, len(diff(invec(C), vA)), len(vA) - 1))
        cur = C
    child = x; last = []
    for kind in 'FR':
        for v in range(1, N + 1):
            a = child if kind == 'F' else h.pathAlgebra.dualPathAlgebra(child)
            if not h.mutation.mutationIsPossibleAtVertex(a, v): continue
            try: ch, J = h.step(child, kind + str(v))
            except AssertionError: continue
            if J or h.search._coxeterKeyOrNone(ch) != pk0: continue
            nm = h.lna_name(ch, N)
            if nm: last.append((kind + str(v), nm, ch))
    print('J = 0 key-kept steps from the J != 0 child that land on an LNA:', [(m, nm) for m, nm, _ in last])
    for mv, nm, ch in last:
        h.check(child, mv, 'last step -> ' + nm)
        print('      LNA %s: invariants == %s: %s; == child: %s' % (nm, B, invec(lnaC(nm)) == vB, invec(lnaC(nm)) == invec(cur)))
    # controls
    print('\ncontrols:')
    ok = sum(invec(mul(mul(R, CA), [list(r) for r in zip(*R)])) == vA for R in (unimod(N, rng) for _ in range(10)))
    print('  YES (code): r C_A r^T, r in GL_12(Z), 10 draws: invariants equal on %d of 10' % ok)
    for dt in (2, 3):
        byinv = byuni = 0
        for _ in range(10):
            R = unimod(N, rng, dt); C2 = mul(mul(R, CA), [list(r) for r in zip(*R)])
            try: byinv += invec(C2) != vA
            except AssertionError: byuni += 1         # det C2 = det^2 != 1: phi() refuses
        print('  planted error: r with det %d, 10 draws: caught by the invariant %d, by the unimodularity check %d' % (dt, byinv, byuni))
    a11, b11 = '003033030', '030330030'
    print('  E-318 n = 11 pair %s / %s separating f != x-1:' % (a11, b11), sorted({k[0] for k in diff(invec(lnaC(a11)), invec(lnaC(b11))) if k != 'charpoly'}))
    print('%.0fs' % (time.time() - t0))
