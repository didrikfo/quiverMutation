"""P/Q separation at n = 10 (round 054): key groups holding >= 2 orbit+mirror classes, and what Hochschild cohomology
and the F-047 Smith profile say about them.  usage: python workshop/rounds/054/maverick_pq.py N [N ...]  [--nosnf]

Hochschild cohomology of an LNA (linear quiver, monomial relations) from the normalised relative bar complex
Hom_{E-E}(r^{(x)m}, A): basis of C^m = chains a_0<...<a_m with consecutive paths nonzero and path (a_0,a_m) nonzero.
Rank over Q by numpy (matrices 0/+-1, at most C(n,m+1) square), d^2 = 0 asserted.
"""
import sys, itertools, collections
import numpy as np
sys.path.insert(0, "workshop/rounds/029")
from quivermutation import coxeterTables as ct, freeMoves

def nonzeroPairs(n, rl):
    """set of (a,b), a<=b, vertices 1..n, whose path is nonzero.  rl[p] = k: relation from p+1 to p+1+k."""
    rels = [(p + 1, p + 1 + k) for p, k in enumerate(rl) if k]
    ok = set()
    for a in range(1, n + 1):
        for b in range(a, n + 1):
            if not any(a <= s and e <= b for s, e in rels):
                ok.add((a, b))
    return ok

def hh(n, rl, ok=None):
    ok = nonzeroPairs(n, rl) if ok is None else ok
    chains = {0: [(a,) for a in range(1, n + 1)]}
    for m in range(1, n):
        L = []
        for c in itertools.combinations(range(1, n + 1), m + 1):
            if all((c[i], c[i + 1]) in ok for i in range(m)) and (c[0], c[-1]) in ok:
                L.append(c)
        if not L: break
        chains[m] = L
    top = max(chains)
    idx = {m: {c: i for i, c in enumerate(chains[m])} for m in chains}
    D = {}
    for m in range(0, top):
        M = np.zeros((len(chains[m + 1]), len(chains[m])))
        for r, c in enumerate(chains[m + 1]):
            # (df)(c) coordinate, f ranges over basis cochains of degree m
            terms = []
            tail = c[1:]
            if tail in idx[m] and (c[0], c[-1]) in ok: terms.append((idx[m][tail], 1))      # x_1 f(a_1..)
            head = c[:-1]
            if head in idx[m] and (c[0], c[-1]) in ok: terms.append((idx[m][head], (-1) ** (m + 1)))
            for i in range(1, m + 1):
                cc = c[:i] + c[i + 1:]
                if cc in idx[m] and (c[0], c[-1]) in ok: terms.append((idx[m][cc], (-1) ** i))
            for j, s in terms: M[r, j] += s
        D[m] = M
    for m in range(0, top - 1):
        assert not np.any(D[m + 1] @ D[m]), "d^2 != 0"
    dims = [len(chains[m]) for m in range(top + 1)]
    rk = [np.linalg.matrix_rank(D[m]) if D[m].size else 0 for m in range(top)]
    out = []
    for m in range(top + 1):
        kern = dims[m] - (rk[m] if m < top else 0)
        im = rk[m - 1] if m > 0 else 0
        out.append(kern - im)
    while out and out[-1] == 0: out.pop()
    return tuple(out)

def classesOf(n):
    keys = ct.lnaKeys(n)
    lnas, orbits = freeMoves.derivedOrbits(n, rules=None, free=True, edges=True, doubles=True)
    oid = {m: o for o, mem in orbits.items() for m in mem}
    par = {o: o for o in orbits}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    closed = collections.defaultdict(bool)
    for m in lnas:
        mr = tuple(freeMoves.mirrorRow(n, m))
        if mr in oid: par[f(oid[m])] = f(oid[mr])
    comps = collections.defaultdict(list)
    for m in lnas: comps[f(oid[m])].append(m)
    # an orbit is mirror-closed when it holds the mirror of its own member
    norb = len(orbits)
    selfm = None
    return keys, comps, oid, norb, selfm

if __name__ == '__main__':
    nosnf = '--nosnf' in sys.argv
    for n in [int(a) for a in sys.argv[1:] if not a.startswith('--')]:
        keys, comps, oid, norb, selfm = classesOf(n)
        hhall = collections.Counter()
        for m in keys: pass
        bykey = collections.defaultdict(list)
        for c, mem in comps.items(): bykey[keys[mem[0]]].append(mem)
        multi = {k: v for k, v in bykey.items() if len(v) > 1}
        print(f"n={n}: LNAs {len(keys)}, move orbits {norb}, orbit+mirror classes {len(comps)}, key groups {len(bykey)}, "
              f"key groups with >=2 classes: {len(multi)}")
        for k, cl in multi.items():
            for mem in cl:
                hhp = collections.Counter(hh(n, m) for m in mem)
                line = f"  key {k[:4]}..  size {len(mem):5d}  orbits {len({oid[m] for m in mem})}  HH profiles {dict(hhp)}"
                if not nosnf:
                    import toolsmith_snfresolve as sr
                    pr = {sr.profile(n, m) for m in mem}
                    line += f"  SNF profiles {len(pr)}"
                    cl_prof = pr
                print(line)
            if not nosnf:
                sets = []
                for mem in cl:
                    import toolsmith_snfresolve as sr
                    sets.append({sr.profile(n, m) for m in mem})
                print("   SNF profile sets equal across classes:", all(s == sets[0] for s in sets), "disjoint:", all(not (s & t) for s in sets for t in sets if s is not t))
