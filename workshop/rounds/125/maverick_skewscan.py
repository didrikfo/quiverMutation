"""Round 125 maverick, companion of maverick_skew.py.  Is X_a = SNF(a C + C^T) = coker(Phi - a) (C unimodular) a function of S = C + C^T?
  plant N    Cartan level: C' = C + K (K integer skew, small, det C' = +-1) has S' = S EXACTLY; count plants with SNF(X_a) changed (a = 3, -3, 2, -2)
             starting from the two E-318 LNAs; also a planted-error control for the congruence test (C -> P^T C Q, Q != P).
  scan N     LNA level: inside every Coxeter-key group at length N, count groups whose SNF(S) classes are split further by
             (SK, X3, X-3, X2, X-2) -- each such split is a pair of same-key, same-SNF(S) LNAs certified derived-inequivalent (H-024 necessity).
  keyplant   Cartan level, Coxeter key KEPT: +-1 changes of 1-2 (all) and 3-4 (RAND random) lower entries of the two E-318 Cartans
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skewscan.py plant 3000 | RAND=20000 ... keyplant | scan 11"""
import sys, os, random, collections, time
sys.path.insert(0, '.')
import flint
from quivermutation import coxeterTables as ct
T0 = time.time()


def fsnf(M):
    n = M.nrows(); a = M.snf(); return tuple(sorted(abs(int(a[i, i])) for i in range(n) if abs(int(a[i, i])) != 1))


def cinv(C):
    T = C.transpose(); out = {'S': fsnf(C + T), 'SK': fsnf(C - T)}
    for a in (3, -3, 2, -2): out['X%d' % a] = fsnf(a * C + T)
    return out


def charpoly(C):
    n = C.nrows(); d = int(C.det()); Ci = flint.fmpz_mat([[int(C.inv()[i, j] * d) for j in range(n)] for i in range(n)]) * d
    return str((-(Ci * C.transpose())).charpoly())


def plant(K):
    rng = random.Random(125); n = 11
    for name in ('003033030', '030330030'):
        C = flint.fmpz_mat(ct.lnaCartanMatrix(n, tuple(int(c) for c in name))); I0 = cinv(C); ch = collections.Counter(); tot = 0
        for _ in range(K):
            Kk = flint.fmpz_mat(n, n)
            for _ in range(rng.randint(1, 3)):
                i, j = rng.sample(range(n), 2); v = rng.choice([-1, 1]); Kk[i, j] += v; Kk[j, i] -= v
            C2 = C + Kk
            if abs(int(C2.det())) != 1: continue
            tot += 1; I2 = cinv(C2); assert I2['S'] == I0['S']
            for k in I0:
                if I2[k] != I0[k]: ch[k] += 1
            if I2['SK'] == I0['SK'] and I2['X3'] != I0['X3']: ch['X3 with SK same'] += 1
            if I2['SK'] == I0['SK'] and I2['X3'] != I0['X3'] and charpoly(C2) == charpoly(C): ch['X3 with SK and Coxeter key same'] += 1
        print('plant %s: %d unimodular C + K (S identical); invariant changed in: %s' % (name, tot, dict(ch)))
    # planted-error control: a non-congruence C -> P^T C Q must be caught by the invariant comparison
    C = flint.fmpz_mat(ct.lnaCartanMatrix(n, (0, 0, 3, 0, 3, 3, 0, 3, 0))); I0 = cinv(C); bad = 0
    def um():
        P = flint.fmpz_mat([[int(i == j) for j in range(n)] for i in range(n)])
        for _ in range(25):
            i, j = rng.sample(range(n), 2); E = flint.fmpz_mat([[int(a == b) for b in range(n)] for a in range(n)]); E[i, j] = rng.choice([-1, 1]); P = P * E
        return P
    for _ in range(100):
        P, Q = um(), um(); bad += cinv(P.transpose() * C * Q) != I0
    print('planted-error control: %d of 100 non-congruent P^T C Q flagged (want > 0)' % bad)


def keyplant(_):
    """key-preserving plants: change 1 or 2 strictly-upper entries of the (unitriangular) LNA Cartan by +-1 (det stays 1); keep those with
    the SAME Coxeter key; tabulate (SNF(S) same?, X3 same?).  A plant with SNF(S) and key equal but X3 different shows X3 is not a function
    of (key, SNF(S)) on unimodular forms."""
    import itertools
    n = 11
    for name in ('003033030', '030330030'):
        C = flint.fmpz_mat(ct.lnaCartanMatrix(n, tuple(int(c) for c in name))); k0 = charpoly(C); I0 = cinv(C)
        low = [(i, j) for i in range(n) for j in range(i)]; assert all(int(C[j, i]) == 0 for i, j in low)   # LNA Cartan is lower unitriangular
        tab = collections.Counter(); ex = None; tot = 0
        rng = random.Random(7); R = int(os.environ.get('RAND', '0'))
        cands = [(pos, sg) for m in (1, 2) for pos in itertools.combinations(low, m) for sg in itertools.product((-1, 1), repeat=m)]
        cands += [(tuple(rng.sample(low, m)), tuple(rng.choice((-1, 1)) for _ in range(m))) for m in (3, 4) for _ in range(R)]
        if True:
            for pos, sg in cands:
                    C2 = flint.fmpz_mat(C)
                    for (i, j), v in zip(pos, sg): C2[i, j] += v
                    if charpoly(C2) != k0: continue
                    tot += 1; I2 = cinv(C2); key = (I2['S'] == I0['S'], I2['X3'] == I0['X3'], I2['SK'] == I0['SK'])
                    tab[key] += 1
                    if key[0] and not key[1] and ex is None: ex = (pos, sg, I0['S'], I0['X3'], I2['X3'])
        print('keyplant %s: %d same-key plants; (S same, X3 same, SK same): %s; example S same / X3 differs: %s  %.0fs'
              % (name, tot, dict(tab), ex, time.time() - T0), flush=True)


def scan(n):
    idx = ct.lnaKeyIndex(n); groups = [g for g in idx.values() if len(g) > 1]
    print('n = %d: %d key groups of size >= 2 (%d LNAs)  %.0fs' % (n, len(groups), sum(map(len, groups)), time.time() - T0), flush=True)
    split = []; nS = nF = 0
    for g in groups:
        I = {r: cinv(flint.fmpz_mat(ct.lnaCartanMatrix(n, tuple(int(c) for c in r)))) for r in g}
        sS = {I[r]['S'] for r in g}; full = {tuple(sorted(I[r].items())) for r in g}; nS += len(sS); nF += len(full)
        if len(full) > len(sS): split.append((len(g), len(sS), len(full), g, I))
    print('S-classes %d, full-invariant classes %d; groups where C-invariants split an SNF(S) class: %d  %.0fs' % (nS, nF, len(split), time.time() - T0))
    for sz, a, b, g, I in split[:12]:
        by = collections.defaultdict(list)
        for r in g: by[I[r]['S']].append(r)
        for s, rs in by.items():
            sub = collections.defaultdict(list)
            for r in rs: sub[tuple((k, v) for k, v in sorted(I[r].items()) if k != 'S')].append(r)
            if len(sub) > 1:
                print('  group %d: SNF(S) %s split into %s by %s; e.g. %s' % (sz, s, [len(v) for v in sub.values()],
                      sorted({k for t in sub for k, v in t if len({dict(u)[k] for u in sub}) > 1}), ' | '.join(ct.className(v[0]) for v in sub.values())))


if __name__ == '__main__':
    {'plant': plant, 'keyplant': keyplant, 'scan': scan}[sys.argv[1]](int(sys.argv[2]) if len(sys.argv) > 2 else 0)
