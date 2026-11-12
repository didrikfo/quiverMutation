"""Round 125 maverick, response to referee: form-level power control for claim (4).
Is X3 = SNF(3C + C^T) a function of (Coxeter key, SNF(S), SNF(SK)) on unimodular lower-unitriangular forms (NOT LNAs)?
Random C, strictly-lower entries in {0, +-1, +-2}, N per n (seed 125).  Counts (key, SNF(S), SNF(SK)) classes with >= 2 X3 values,
and how many of those differ 2-adically (the 2-parts of the X3 tuples differ).  Also checks the referee's n = 4 pair, SNF cross-checked by pari.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_formctl.py [N=60000]"""
import sys, random, collections, time
import flint, cypari
pari = cypari.pari; T0 = time.time()


def fsnf(M):
    a = M.snf(); return tuple(sorted(abs(int(a[i, i])) for i in range(M.nrows()) if abs(int(a[i, i])) != 1))


def v2part(t):
    out = []
    for x in t:
        p = 1
        while x and x % 2 == 0: x //= 2; p *= 2
        out.append(0 if x == 0 else p)
    return tuple(sorted(q for q in out if q != 1))


def inv(C):
    T = C.transpose(); Ci = C.inv()                                   # C unitriangular, inverse integral
    Ci = flint.fmpz_mat([[int(Ci[i, j]) for j in range(C.nrows())] for i in range(C.nrows())])
    return str((-(Ci * T)).charpoly()), fsnf(C + T), fsnf(C - T), fsnf(3 * C + T)


def pari_snf(M):
    n = M.nrows(); d = pari.matsnf(pari.matrix(n, n, [int(M[i, j]) for i in range(n) for j in range(n)]))
    return tuple(sorted(abs(int(x)) for x in d if abs(int(x)) != 1))


def main(N):
    C1 = flint.fmpz_mat([[1, 0, 0, 0], [-2, 1, 0, 0], [1, -2, 1, 0], [-1, 0, 2, 1]])
    C2 = flint.fmpz_mat([[1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 1, 0], [1, -1, 0, 1]])
    for C in (C1, C2):
        k, s, sk, x3 = inv(C); assert x3 == pari_snf(3 * C + C.transpose()) and s == pari_snf(C + C.transpose())
        print('referee pair: key %s  SNF(S) %s  SNF(SK) %s  X3 %s (pari agrees)' % (k, s, sk, x3))
    rng = random.Random(125)
    for n in (4, 5, 6):
        cl = collections.defaultdict(dict)
        for _ in range(N):
            C = flint.fmpz_mat([[1 if i == j else (rng.choice((0, 0, 1, -1, 2, -2)) if j < i else 0) for j in range(n)] for i in range(n)])
            k, s, sk, x3 = inv(C); cl[(k, s, sk)].setdefault(x3, C)
        hit = {c: d for c, d in cl.items() if len(d) > 1}
        two = {c: d for c, d in hit.items() if len({v2part(x) for x in d}) > 1}
        print('n = %d: %d classes (key, SNF(S), SNF(SK)); %d with >= 2 values of X3; %d of them differ 2-adically  %.0fs'
              % (n, len(cl), len(hit), len(two), time.time() - T0), flush=True)
        if two:
            c, d = min(two.items(), key=lambda t: len(t[0][0]))
            print('   2-adic example: key %s SNF(S) %s SNF(SK) %s; X3 values %s' % (c[0], c[1], c[2], sorted(d)))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 60000)
