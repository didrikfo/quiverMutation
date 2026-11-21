"""Round 134 scholar: object-level nu^m S_l on simple modules of a linear Nakayama algebra (E-339 blocker, H-025).
S_l is given by its minimal projective resolution; nu is applied by E-246's minimal-bimodule method
(workshop/rounds/089/maverick_bimod.py, unchanged) and the complex is minimised after every step.
For X = nu^m S_l we print: number of indecomposable projective terms (rank), total cohomology dimension
sum_i dim H^*(e_i X) (computed as H^*(Hom(P_i, X))), the cohomological degree range, and the Euler vector.
Exact criterion: nu^m S_l ~ S_j[k] for some j, k  iff  total cohomology dimension = 1 (a complex with cohomology a
simple module in one degree is isomorphic to it in D^b).  FCY (nu^a ~ [b] as functors, or just on S_l) forces the
total cohomology dim of nu^m S_l to be bounded in m; linear growth along m = order, 2 order, 3 order is evidence
against (exact at each m, the boundedness itself is not proved by a finite sample).
Consistency check: Euler(nu^m S_l) must equal M^m e_l, M = Dmat C^{-1} the matrix of nu on K0 from the Cartan data alone
(independent of the bimodule resolution Q).  --plant deletes one coefficient of Q's differential (planted error).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/134/scholar_simple_nu.py RELSTRING M1,M2,... [l1,l2,...] [--plant]"""
import sys, os, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "089")); sys.path.insert(0, os.path.join(HERE, "..", "087"))
from maverick_bimod import resolve, nu_bimod, rank_mod
from scholar_rankgrowth import Alg, Cx, P

def simple_resolution(A, l):
    """minimal projective resolution of S_l: P_{x_k} in degree -k, canonical maps (coefficient 1)"""
    # kernel of P_{x_k} -> P_{x_{k-1}} has composition factors x_{k-1} + p_{x_{k-1}} .. x_k + p_{x_k} - 1
    xs = [l]
    if A.p[l] >= 2:
        xs.append(l + 1)
        while True:
            a, b = xs[-2], xs[-1]
            st, sp = a + A.p[a], b + A.p[b] - 1
            if st > sp: break
            xs.append(st)
    X = Cx()
    for k, x in enumerate(xs): X.add(k, -k, x)
    for k in range(1, len(xs)): X.setd(k, k - 1, 1)
    return X

def cohomology(A, X):
    """{(vertex i, degree t): dim H^t(Hom(P_i, X))}, nonzero entries only; asserts d^2 = 0"""
    out = {}
    for i in range(1, A.n + 1):
        vs = [u for u, (t, a) in X.vert.items() if A.hom(i, a)]
        if not vs: continue
        S = set(vs); deg = {u: X.vert[u][0] for u in vs}
        ent = {(u, v): c for u in vs for v, c in X.out.get(u, {}).items() if v in S}
        bys = collections.defaultdict(dict)
        for (u, v), c in ent.items(): bys[u][v] = c
        for u in bys:
            acc = collections.defaultdict(int)
            for v, c in bys[u].items():
                for w, c2 in bys.get(v, {}).items():
                    if A.hom(i, X.vert[w][1]): acc[w] = (acc[w] + c * c2) % P
            assert all(x == 0 for x in acc.values()), "d^2 != 0"
        degs = sorted(set(deg.values())); rk = {}
        for dg in degs:
            src = [u for u in vs if deg[u] == dg]; tgt = [u for u in vs if deg[u] == dg + 1]
            rk[dg] = rank_mod([[ent.get((s, t), 0) for s in src] for t in tgt], len(src)) if src and tgt else 0
        for dg in degs:
            h = sum(1 for u in vs if deg[u] == dg) - rk[dg] - rk.get(dg - 1, 0)
            if h: out[(i, dg)] = h
    return out

def euler(A, H):
    v = [0] * (A.n + 1)
    for (i, t), h in H.items(): v[i] += (-1) ** (t % 2) * h
    return v[1:]

def main():
    rel = sys.argv[1]; ms = sorted(int(x) for x in sys.argv[2].split(","))
    A = Alg([int(c) for c in rel]); n = A.n
    pos = [a for a in sys.argv[3:] if not a.startswith("--")]
    ls = [int(x) for x in pos[0].split(",")] if pos else list(range(1, n + 1))
    t0 = time.time(); levels = resolve(A)
    if "--plant" in sys.argv:     # planted error: delete the first nonzero coefficient of Q's differential
        for lv in levels[1:]:
            d = next((d for d in lv["diff"] if d), None)
            if d: g = next(iter(d)); del d[g]; print("PLANTED: one coefficient of Q differential deleted"); break
    # nu on K0, independent of Q: [nu P_j] = [I_j], so M = Dmat C^{-1} (C: columns dim P_j; Dmat: columns dim I_j)
    from fractions import Fraction
    C = [[Fraction(int(A.hom(i, j))) for j in range(1, n + 1)] for i in range(1, n + 1)]   # C[i][j] = [P_j : S_i]
    Dm = [[Fraction(int(A.hom(j, i))) for j in range(1, n + 1)] for i in range(1, n + 1)]  # [I_j : S_i] = dim Hom(P_i,I_j)... = hom(j,i)
    aug = [C[i][:] + [Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        pr = next(r for r in range(c, n) if aug[r][c]); aug[c], aug[pr] = aug[pr], aug[c]
        aug[c] = [x / aug[c][c] for x in aug[c]]
        for r in range(n):
            if r != c and aug[r][c]: f = aug[r][c]; aug[r] = [x - f * y for x, y in zip(aug[r], aug[c])]
    Ci = [row[n:] for row in aug]
    M = [[sum(Dm[i][k] * Ci[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    assert all(x.denominator == 1 for r in M for x in r)
    M = [[int(x) for x in r] for r in M]
    def matvec(v): return [sum(M[i][j] * v[j] for j in range(n)) for i in range(n)]
    print(f"{rel} n={n} m={ms} simples={ls}  (Q built {time.time()-t0:.1f}s)", flush=True)
    for l in ls:
        X = simple_resolution(A, l); H0 = cohomology(A, X)
        assert H0 == {(l, 0): 1}, ("bad resolution", l, H0)
        e = [0] * n; e[l - 1] = 1; k = 0
        for m in ms:
            while k < m:
                X, _, viol = nu_bimod(A, levels, X); e = matvec(e); k += 1
            H = cohomology(A, X); eu = euler(A, H); tot = sum(H.values())
            degs = [t for (_, t) in H]
            shift = None
            if tot == 1: (j, t), = H.keys(); shift = f"S_{j}[{-t}]"
            print(f"l={l:2d} m={m:4d} rank={len(X.vert):5d} H-total={tot:5d} deg=[{min(degs)},{max(degs)}] "
                  f"Euler=M^m e_l:{eu == e} Euler={'+e_l' if eu == e_l(n, l, 1) else '-e_l' if eu == e_l(n, l, -1) else 'other'} "
                  f"shift-of-simple:{shift}  {time.time()-t0:.0f}s", flush=True)

def e_l(n, l, s):
    v = [0] * n; v[l - 1] = s; return v

if __name__ == "__main__":
    main()
