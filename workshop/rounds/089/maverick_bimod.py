"""Round 089 maverick: rank growth of nu^k P_l via the MINIMAL projective bimodule resolution Q of D(Lambda).
nu X = X (x)_Lambda Q (totalised): strict d^2 = 0, no bar resolution, no homotopy transfer; then Gaussian minimisation
(minimise() of scholar_rankgrowth).  Q is computed by linear algebra over GF(p) on the category C (objects 1..n,
C(a,b) one-dim iff hom(a,b)).  usage (repo root): timeout 10m .venv/bin/python workshop/rounds/089/maverick_bimod.py RELSTRING K [--check]"""
import sys, os, collections, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "087"))
from scholar_rankgrowth import Alg, Cx, minimise, P

def solve_cols(basis, vec):
    """coordinates of vec in the list of independent vectors basis (None if not in span)"""
    m = len(basis)
    if m == 0: return [] if not any(vec) else None
    L = len(vec)
    rows = [[basis[j][i] for j in range(m)] + [vec[i]] for i in range(L)]
    piv = []; r = 0
    for c in range(m):
        pr = next((i for i in range(r, L) if rows[i][c] % P), None)
        if pr is None: return None
        rows[r], rows[pr] = rows[pr], rows[r]
        iv = pow(rows[r][c], P - 2, P); rows[r] = [x * iv % P for x in rows[r]]
        for i in range(L):
            if i != r and rows[i][c] % P:
                f = rows[i][c]; rows[i] = [(x - f * y) % P for x, y in zip(rows[i], rows[r])]
        r += 1
    if any(rows[i][m] % P for i in range(r, L)): return None
    return [rows[i][m] for i in range(m)]
def nullspace(M, ncols):
    """M list of rows; basis of {x : M x = 0}"""
    rows = [r[:] for r in M]; piv = []; r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(rows)) if rows[i][c] % P), None)
        if pr is None: continue
        rows[r], rows[pr] = rows[pr], rows[r]
        iv = pow(rows[r][c], P - 2, P); rows[r] = [x * iv % P for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] % P:
                f = rows[i][c]; rows[i] = [(x - f * y) % P for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    out = []
    for fc in range(ncols):
        if fc in piv: continue
        v = [0] * ncols; v[fc] = 1
        for i, pc in enumerate(piv): v[pc] = (-rows[i][fc]) % P
        out.append(v)
    return out
def span_reduce(vecs):
    """independent subset basis (echelon list) of vecs"""
    ech = []; 
    for v in vecs:
        w = v[:]
        for (pc, e) in ech:
            if w[pc] % P:
                f = w[pc]; w = [(x - f * y) % P for x, y in zip(w, e)]
        pc = next((i for i, x in enumerate(w) if x % P), None)
        if pc is None: continue
        iv = pow(w[pc], P - 2, P); w = [x * iv % P for x in w]; ech.append((pc, w))
    return ech
def matvec(M, v): return [sum(a * b for a, b in zip(r, v)) % P for r in M]

class BiMod:
    """dim[(a,x)], A[(a,x)] : M(a,x)->M(a-1,x) (rows=target), X[(a,x)] : M(a,x)->M(a,x+1)"""
    def __init__(s, dim, A, X): s.dim, s.A, s.X = dim, A, X

def resolve(A):
    n = A.n; hom = A.hom
    verts = [(a, x) for a in range(1, n + 1) for x in range(1, n + 1)]
    dim = {v: (1 if hom(*v) else 0) for v in verts}
    ma = {}; mx = {}
    for (a, x) in verts:
        if a > 1 and hom(a, a - 1):
            ma[(a, x)] = [[1 if (dim[(a, x)] and dim[(a - 1, x)]) else 0] * dim[(a, x)] for _ in range(dim[(a - 1, x)])]
        if x < n and hom(x + 1, x):
            mx[(a, x)] = [[1 if (dim[(a, x)] and dim[(a, x + 1)]) else 0] * dim[(a, x)] for _ in range(dim[(a, x + 1)])]
    M = BiMod(dim, ma, mx)
    levels = []      # per level: list of gens (c,d, vec_in_prev_gens_basis) ; vec: dict prevgen index -> coeff
    prevconv = None
    while sum(M.dim.values()) > 0:
        # top
        gens = []
        for (a, x) in verts:
            if M.dim[(a, x)] == 0: continue
            img = []
            if a + 1 <= n and (a + 1, x) in M.A:
                T = M.A[(a + 1, x)]; img += [[T[i][j] for i in range(len(T))] for j in range(M.dim[(a + 1, x)])]
            if x - 1 >= 1 and (a, x - 1) in M.X:
                T = M.X[(a, x - 1)]; img += [[T[i][j] for i in range(len(T))] for j in range(M.dim[(a, x - 1)])]
            ech = span_reduce(img)
            for k in range(M.dim[(a, x)]):
                e = [0] * M.dim[(a, x)]; e[k] = 1
                ech2 = span_reduce([w for _, w in ech] + [e])
                if len(ech2) > len(ech): ech = ech2; gens.append((a, x, e))
        # free module F: basis at (a,x) = gens g with hom(c,a), hom(x,d)
        def fbasis(a, x): return [gi for gi, (c, d, _) in enumerate(gens) if hom(c, a) and hom(x, d)]
        FB = {v: fbasis(*v) for v in verts}
        def act_gen(gi, a, x):
            c, d, w = gens[gi]; vec = w[:]; cur = c
            while cur > a: vec = matvec(M.A[(cur, d)], vec) if (cur, d) in M.A else [0] * M.dim[(cur - 1, d)]; cur -= 1
            cd = d
            while cd < x: vec = matvec(M.X[(a, cd)], vec) if (a, cd) in M.X else [0] * M.dim[(a, cd + 1)]; cd += 1
            return vec
        # pi and kernel
        Kb = {}
        for v in verts:
            a, x = v; fb = FB[v]
            if not fb: Kb[v] = []; continue
            cols = [act_gen(gi, a, x) for gi in fb]
            mat = [[cols[j][i] for j in range(len(fb))] for i in range(M.dim[v])]
            Kb[v] = nullspace(mat, len(fb)) if M.dim[v] else [[1 if i == j else 0 for i in range(len(fb))] for j in range(len(fb))]
        diff = []
        for (c, d, w) in gens:
            if prevconv is None: diff.append({}); continue
            pFB, pKb = prevconv; vec = [0] * len(pFB[(c, d)])
            for wi, kb in zip(w, pKb[(c, d)]):
                if wi:
                    for t, kv in enumerate(kb): vec[t] = (vec[t] + wi * kv) % P
            diff.append({pFB[(c, d)][t]: vec[t] for t in range(len(vec)) if vec[t]})   # prev gen index -> coeff
        levels.append({"gens": [(c, d) for c, d, _ in gens], "diff": diff})
        # new module K
        kdim = {v: len(Kb[v]) for v in verts}
        # F actions: basis idx lists
        def f_a(a, x, vec):    # F(a,x) -> F(a-1,x), vec in basis FB[(a,x)]
            src = FB[(a, x)]; tgt = FB[(a - 1, x)]; out = [0] * len(tgt)
            for coef, gi in zip(vec, src):
                if gi in tgt: out[tgt.index(gi)] = (out[tgt.index(gi)] + coef) % P
            return out
        def f_x(a, x, vec):
            src = FB[(a, x)]; tgt = FB[(a, x + 1)]; out = [0] * len(tgt)
            for coef, gi in zip(vec, src):
                if gi in tgt: out[tgt.index(gi)] = (out[tgt.index(gi)] + coef) % P
            return out
        nA = {}; nX = {}
        for (a, x) in verts:
            if kdim[(a, x)] == 0: continue
            if a > 1 and kdim[(a - 1, x)] and hom(a, a - 1):
                cols = []
                for b in Kb[(a, x)]:
                    co = solve_cols(Kb[(a - 1, x)], f_a(a, x, b)); assert co is not None; cols.append(co)
                nA[(a, x)] = [[cols[j][i] for j in range(len(cols))] for i in range(kdim[(a - 1, x)])]
            if x < n and kdim[(a, x + 1)] and hom(x + 1, x):
                cols = []
                for b in Kb[(a, x)]:
                    co = solve_cols(Kb[(a, x + 1)], f_x(a, x, b)); assert co is not None; cols.append(co)
                nX[(a, x)] = [[cols[j][i] for j in range(len(cols))] for i in range(kdim[(a, x + 1)])]
        prevconv = (FB, Kb)
        M = BiMod(kdim, nA, nX)
    return levels

def nu_bimod(A, levels, X):
    """X : Cx of projectives; returns minimised Cx for nu X"""
    hom = A.hom; T = Cx(); idx = {}
    for u, (t, a) in X.vert.items():
        for j, lv in enumerate(levels):
            for gi, (c, d) in enumerate(lv["gens"]):
                if hom(c, a): idx[(u, j, gi)] = len(idx); T.add(idx[(u, j, gi)], t - j, d)
    for (u, j, gi), k in idx.items():
        t, a = X.vert[u]
        # Q-differential: gen gi of level j maps to prev gens gp with coeff kappa: P_d -> P_{d'}, needs hom(c',a)
        if j >= 1:
            for gp, kap in levels[j]["diff"][gi].items():
                tgt = idx.get((u, j - 1, gp))
                if tgt is not None: T.setd(k, tgt, kap)
        for v, cf in X.out.get(u, {}).items():
            tgt = idx.get((v, j, gi))
            if tgt is not None: T.setd(k, tgt, cf * (-1) ** j)
    Y, viol = minimise(A, T)
    Z = Cx(); ren = {}
    for u, dv in Y.vert.items(): ren[u] = len(ren); Z.vert[ren[u]] = dv
    for u in Y.vert:
        for v, c in Y.out[u].items(): Z.out[ren[u]][ren[v]] = c; Z.inn[ren[v]][ren[u]] = c
    return Z, len(T.vert), viol

def run(rl, K, verbose=True):
    A = Alg(rl); t0 = time.time(); levels = resolve(A)
    if verbose: print("n", A.n, "p", A.p[1:], "Q gens per level", [len(l["gens"]) for l in levels], f"{time.time()-t0:.1f}s", flush=True)
    cur = {}
    for l in range(1, A.n + 1):
        X = Cx(); X.vert[0] = (0, l); cur[l] = X
    res = [[1] * A.n]
    for k in range(1, K + 1):
        nxt = {}; big = 0
        for l, X in cur.items():
            Z, size, v = nu_bimod(A, levels, X); nxt[l] = Z; big = max(big, size)
        cur = nxt; ranks = [len(X.vert) for X in cur.values()]; res.append(ranks)
        if verbose: print(f"k={k}: R_k {sum(ranks)} max {max(ranks)} biggest {big} {time.time()-t0:.0f}s", flush=True)
    return A, res, cur
if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    A, res, cur = run([int(c) for c in args[0]], int(args[1]))
    print("R_k", [sum(r) for r in res])

def rank_mod(rows, ncols):
    rows = [r[:] for r in rows]; r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(rows)) if rows[i][c] % P), None)
        if pr is None: continue
        rows[r], rows[pr] = rows[pr], rows[r]; iv = pow(rows[r][c], P - 2, P); rows[r] = [x * iv % P for x in rows[r]]
        for i in range(r + 1, len(rows)):
            if rows[i][c] % P:
                f = rows[i][c]; rows[i] = [(x - f * y) % P for x, y in zip(rows[i], rows[r])]
        r += 1
    return r
def hom_total(A, X, i, into):
    """total dim H^*(Hom(P_i, X)) (into=False) or H^*(Hom(X, P_i)) (into=True); also checks d^2 = 0 on the Hom complex"""
    hom = A.hom
    vs = [u for u, (t, a) in X.vert.items() if (hom(a, i) if into else hom(i, a))]
    pos = {u: k for k, u in enumerate(vs)}
    deg = {u: (-X.vert[u][0] if into else X.vert[u][0]) for u in vs}
    # differential matrix entries (src->tgt)
    ent = {}
    for u in vs:
        for v, c in X.out.get(u, {}).items():
            if v in pos:
                if into: ent[(pos[v], pos[u])] = c     # transpose: v -> u
                else: ent[(pos[u], pos[v])] = c
    # d^2 = 0 check
    byS = collections.defaultdict(dict)
    for (s, t), c in ent.items(): byS[s][t] = c
    for s in byS:
        acc = collections.defaultdict(int)
        for t, c in byS[s].items():
            for w, c2 in byS.get(t, {}).items(): acc[w] = (acc[w] + c * c2) % P
        assert all(x == 0 for x in acc.values()), "d^2 != 0 on Hom complex"
    N = len(vs); degs = sorted(set(deg.values())); tot = 0
    rk = {}
    for dg in degs:
        src = [u for u in vs if deg[u] == dg]; tgt = [u for u in vs if deg[u] == dg + 1]
        mat = [[ent.get((pos[s], pos[t]), 0) for s in src] for t in tgt]
        rk[dg] = rank_mod(mat, len(src)) if src and tgt else 0
    for dg in degs:
        n_dg = sum(1 for u in vs if deg[u] == dg)
        tot += n_dg - rk[dg] - rk.get(dg - 1, 0)
    return tot
