"""Round 135 maverick: a second implementation of nu (Serre functor) on D^b(mod A), A a linear Nakayama algebra.

Method (independent of rounds/089 maverick_bimod and rounds/087 scholar_rankgrowth; shares no code with them):
  X a bounded complex of projectives (each P_i = interval [i, U_i]); Hom(P_i, P_j) is <= 1-dim, so X is
  "types + scalar matrices". nu is applied TERMWISE: P_i -> I_i = [L_i, i] and the standard map P_i -> P_j
  (identity on the overlap) -> the standard map I_i -> I_j. (On the additive category of projectives this IS nu
  up to rescaling the isomorphisms nu(P_i) = I_i: every Hom is spanned by a composite of arrows, so the
  multiplicative scalars can be absorbed.) The complex of injectives is then replaced by a quasi-isomorphic
  complex of projectives by the cone algorithm, at module level (representations, vertexwise linear algebra):
  descending in degree, W^k = cone cycles = {(x, y) in P^{k+1} + J^k : d x = 0, f x = d y}, P^k := projective
  cover of W^k, d := -x-part, f := y-part.  Then Gaussian elimination of iso entries makes the complex minimal.
  No bimodule resolution, no bar resolution.
Checks: every map asserted to be a module map; d^2 = 0; H(P) = H(J) vertexwise (quasi-iso dims); optional
  Serre-duality check dim H^t Hom(X, P_i) = dim H^{-t}(nu X)_i (--serre).
Wrong nu's (controls): --plant=scale  scales every standard injective map along a path of length 2 by 2 (not a
  functor; d^2 = 0 survives);  --plant=kup  uses I_i = [i - p_i + 1, i] (projective length mistaken for injective
  length: a functor to modules, wrong objects);  --plant=trunc  I_i := [L_i + 1, i] (one vertex short);  --plant=long  I_i := [L_i - 1, i];
  --plant=id  nu := identity (a genuine autoequivalence, not Serre).
Usage (repo root): timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py REL m1,m2,.. [objs] [--prime p] [--plant=scale|kup|trunc|long|id] [--serre]
  REL = relation-length digits, n = len(REL) + 2 (same naming as the record). objs: S (simples), P (projectives),
  U (all uniserials), or e.g. S5,P2,U4-5 ; default S,P.
"""
import sys, time
import numpy as np

P = 2147483647
PLANT = False

def args():
    global P, PLANT
    a = sys.argv[1:]
    if "--prime" in a:
        i = a.index("--prime"); P = int(a[i + 1]); del a[i:i + 2]
    for md in ("scale", "kup", "id", "trunc", "long"):
        if "--plant=" + md in a: PLANT = md; a.remove("--plant=" + md)
    serre = "--serre" in a
    if serre: a.remove("--serre")
    return a, serre

# ---------------- linear algebra mod P (int64; P < 2^31 so products fit) ----------------
def rref(M):
    M = np.array(M, dtype=np.int64) % P
    r, c = M.shape; piv = []; row = 0
    for col in range(c):
        if row >= r: break
        nz = np.nonzero(M[row:, col])[0]
        if len(nz) == 0: continue
        k = row + nz[0]
        if k != row: M[[row, k]] = M[[k, row]]
        inv = pow(int(M[row, col]), P - 2, P)
        M[row] = (M[row] * inv) % P
        others = np.nonzero(M[:, col])[0]
        for o in others:
            if o != row: M[o] = (M[o] - M[o, col] * M[row]) % P
        piv.append(col); row += 1
    return M[:row], piv

def mulmod(A, B):
    """A @ B mod P without int64 overflow (split A into 16-bit halves)"""
    A = np.asarray(A, dtype=np.int64) % P; B = np.asarray(B, dtype=np.int64) % P
    hi, lo = A >> 16, A & 0xFFFF
    return (((hi @ B) % P) * 65536 + (lo @ B)) % P

def rank(M):
    M = np.asarray(M)
    if M.size == 0: return 0
    return len(rref(M)[1])

def nullspace(M, ncols):
    """columns spanning ker M (M: rows x ncols)"""
    if ncols == 0: return np.zeros((0, 0), dtype=np.int64)
    if M.shape[0] == 0: return np.eye(ncols, dtype=np.int64)
    R, piv = rref(M); free = [j for j in range(ncols) if j not in piv]
    N = np.zeros((ncols, len(free)), dtype=np.int64)
    for t, f in enumerate(free):
        N[f, t] = 1
        for i, pc in enumerate(piv): N[pc, t] = (-R[i, f]) % P
    return N

def extend(B, W):
    """columns of W (spanning a space containing span B) that extend B to a basis of span W"""
    if W.shape[1] == 0: return []
    cur = B.shape[1] if B.size else 0; rk = rank(B.T) if cur else 0; out = []
    S = B.copy() if cur else np.zeros((W.shape[0], 0), dtype=np.int64)
    for j in range(W.shape[1]):
        T = np.concatenate([S, W[:, j:j + 1]], axis=1)
        r2 = rank(T.T)
        if r2 > rk: S = T; rk = r2; out.append(W[:, j])
    return out

# ---------------- the algebra ----------------
class Alg:
    def __init__(self, rel):
        r = [int(c) for c in rel]; n = len(r) + 2; self.n = n
        # U_i = last vertex of P_i: paths from i stop before a zero relation j -> j + r_j contained in them
        self.U = [None] + [min([n] + [j + r[j - 1] - 1 for j in range(i, n + 1) if j - 1 < len(r) and r[j - 1] > 0]) for i in range(1, n + 1)]
        for i in range(1, n): assert self.U[i] <= self.U[i + 1]
        self.L = [None] + [min(s for s in range(1, a + 1) if self.U[s] >= a) for a in range(1, n + 1)]
    def proj(self, i): return (i, self.U[i])
    def inj(self, i): return (self.L[i], i)

def std_ok(s, t):
    """standard map interval s=(a,b) -> t=(c,d) (identity on overlap) is a nonzero module map iff c <= a <= d <= b
    (arrows v -> v+1; image of s is a quotient [a, d] of s and a submodule of t)."""
    return t[0] <= s[0] <= t[1] <= s[1]

# a complex: deg -> list of interval types; D[k] : matrix rows = deg k+1 summands, cols = deg k summands
class Cx:
    def __init__(self): self.T = {}; self.D = {}
    def degs(self): return sorted(k for k in self.T if self.T[k])
    def mat(self, k):
        r, c = len(self.T.get(k + 1, [])), len(self.T.get(k, []))
        M = self.D.get(k)
        if M is None or M.shape != (r, c): M = np.zeros((r, c), dtype=np.int64); self.D[k] = M
        return M

def vbasis(types, v): return [a for a, t in enumerate(types) if t[0] <= v <= t[1]]

def vmat(X, k, v):
    """differential d^k at vertex v as a matrix (identity on overlaps)"""
    src = vbasis(X.T.get(k, []), v); tgt = vbasis(X.T.get(k + 1, []), v); M = X.mat(k)
    return M[np.ix_(tgt, src)] if src and tgt else np.zeros((len(tgt), len(src)), dtype=np.int64)

def check(X):
    for k in X.degs():
        M = X.mat(k)
        for b, a in zip(*np.nonzero(M)):
            assert std_ok(X.T[k][a], X.T[k + 1][b]), ("not a module map", X.T[k][a], X.T[k + 1][b])
        if k + 1 in X.T and X.T[k + 1]:
            for v in range(1, NV + 1):
                A1 = vmat(X, k, v); A2 = vmat(X, k + 1, v)
                if A1.size and A2.size: assert not mulmod(A2, A1).any(), ("d^2 != 0", k, v)

def cohom(X):
    """{(v, t): dim H^t(X)_v} = dims of H^t Hom(P_v, X)"""
    H = {}
    ks = X.degs()
    if not ks: return H
    for v in range(1, NV + 1):
        for t in range(ks[0], ks[-1] + 1):
            dim = len(vbasis(X.T.get(t, []), v))
            if not dim: continue
            h = dim - rank(vmat(X, t, v)) - rank(vmat(X, t - 1, v))
            if h: H[(v, t)] = h
    return H

def minimise(X):
    changed = True
    while changed:
        changed = False
        for k in X.degs():
            M = X.mat(k); Tk, Tk1 = X.T[k], X.T.get(k + 1, [])
            hit = None
            for b, a in zip(*np.nonzero(M)):
                if Tk[a] == Tk1[b]: hit = (b, a); break
            if hit is None: continue
            b, a = hit; inv = pow(int(M[b, a]), P - 2, P)
            col = M[:, a].copy(); rowv = M[b, :].copy()
            upd = (np.outer(col, rowv) % P) * inv % P
            # composite c -> b ~ a -> r survives only where the standard map c -> r is a module map
            lo0 = np.array([t[0] for t in Tk]); hi0 = np.array([t[1] for t in Tk])
            lo1 = np.array([t[0] for t in Tk1]); hi1 = np.array([t[1] for t in Tk1])
            ok = (lo1[:, None] <= lo0[None, :]) & (lo0[None, :] <= hi1[:, None]) & (hi1[:, None] <= hi0[None, :])
            upd = upd * ok
            M = (M - upd) % P
            keepc = [c for c in range(len(Tk)) if c != a]; keepr = [r for r in range(len(Tk1)) if r != b]
            X.D[k] = M[np.ix_(keepr, keepc)]
            if k - 1 in X.T: X.D[k - 1] = X.mat(k - 1)[keepc, :]
            if k + 1 in X.D and k + 2 in X.T: X.D[k + 1] = X.mat(k + 1)[:, keepr]
            X.T[k] = [Tk[c] for c in keepc]; X.T[k + 1] = [Tk1[r] for r in keepr]
            changed = True; break
    return X

# ---------------- projective replacement (cone algorithm) ----------------
def arrow(types, v):
    """arrow v -> v+1 on a direct sum of intervals: matrix (basis at v+1) x (basis at v)"""
    s = vbasis(types, v); t = vbasis(types, v + 1); A = np.zeros((len(t), len(s)), dtype=np.int64)
    for i, a in enumerate(s):
        if a in t: A[t.index(a), i] = 1
    return A

def replace(J, alg):
    """J: Cx of interval modules (std maps). Returns Cx of projectives quasi-isomorphic to J."""
    n = NV; ks = J.degs(); b = ks[-1]; a = ks[0]
    X = Cx(); f = {}            # f[k][v] : matrix J^k_v x P^k_v
    k = b
    while True:
        Pk1 = X.T.get(k + 1, []); Jk = J.T.get(k, []); Jk1 = J.T.get(k + 1, [])
        newT = []; xcols = []; ycols = []      # generator data
        Wprev = None
        for v in range(1, n + 1):
            pb = vbasis(Pk1, v); jb = vbasis(Jk, v)
            E = len(pb) + len(jb)
            if E == 0: Wprev = None; continue
            # phi(x, y) = (d_P x, f x - d_J y) into P^{k+2}_v + J^{k+1}_v
            dP = vmat(X, k + 1, v) if (k + 2) in X.T else np.zeros((0, len(pb)), dtype=np.int64)
            fv = f.get(k + 1, {}).get(v, np.zeros((len(vbasis(Jk1, v)), len(pb)), dtype=np.int64))
            dJ = vmat(J, k, v)
            top = np.concatenate([dP, np.zeros((dP.shape[0], len(jb)), dtype=np.int64)], axis=1)
            bot = np.concatenate([fv, (-dJ) % P], axis=1)
            Phi = np.concatenate([top, bot], axis=0) % P
            W = nullspace(Phi, E)
            # image of W_{v-1} under the arrow of E
            if Wprev is not None and Wprev.shape[1]:
                Ar = np.zeros((E, Wprev.shape[0]), dtype=np.int64)
                pbm = vbasis(Pk1, v - 1); jbm = vbasis(Jk, v - 1)
                Ap = arrow(Pk1, v - 1); Aj = arrow(Jk, v - 1)
                Ar[:len(pb), :len(pbm)] = Ap; Ar[len(pb):, len(pbm):] = Aj
                B = (Ar @ Wprev) % P
            else: B = np.zeros((E, 0), dtype=np.int64)
            for w in extend(B, W):
                newT.append(alg.proj(v)); xcols.append((v, pb, w[:len(pb)])); ycols.append((v, jb, w[len(pb):]))
            Wprev = W
        if not newT and k < a: break
        X.T[k] = newT
        if Pk1:
            D = np.zeros((len(Pk1), len(newT)), dtype=np.int64)
            for c, (v, pb, x) in enumerate(xcols):
                for i, g in enumerate(pb): D[g, c] = (-x[i]) % P
            X.D[k] = D
        # f^k at each vertex: generator y propagated by arrows of J^k
        fk = {}
        for v in range(1, n + 1):
            pv = vbasis(newT, v); jv = vbasis(Jk, v); fk[v] = np.zeros((len(jv), len(pv)), dtype=np.int64)
        for c, (v0, jb, y) in enumerate(ycols):
            vec = y.copy(); v = v0
            while True:
                pv = vbasis(newT, v)
                if c not in pv: break
                if len(vec): fk[v][:, pv.index(c)] = vec
                if v + 1 > n: break
                vec = (arrow(Jk, v) @ vec) % P
                v += 1
        f[k] = fk
        if not newT and k < a: break
        k -= 1
        assert k > a - 3 * n - 5, "no termination"
    for kk in list(X.T):
        if not X.T[kk]: del X.T[kk]
    return X

def inj(i):
    if PLANT == "kup": return (max(1, i - ALG.U[i] + i), i)     # bug: injective length taken = projective length
    if PLANT == "trunc": L, _ = ALG.inj(i); return (min(L + 1, i), i)   # bug: injective envelope missing its top vertex
    if PLANT == "long": L, _ = ALG.inj(i); return (max(1, L - 1), i)      # bug: injective one vertex too long
    return ALG.inj(i)

def nu_termwise(X):
    J = Cx()
    if PLANT == "id": J.T = {k: list(ts) for k, ts in X.T.items()}; J.D = {k: X.mat(k).copy() for k in X.degs()}; return J
    for k, ts in X.T.items(): J.T[k] = [inj(t[0]) for t in ts]
    for k in X.degs():
        if k + 1 not in X.T: continue
        M = X.mat(k).copy()
        if PLANT == "scale":
            for b, a in zip(*np.nonzero(M)):
                if X.T[k][a][0] - X.T[k + 1][b][0] == 2: M[b, a] = 2 * M[b, a] % P
        J.D[k] = M
    check(J)
    return J

def nu(X):
    J = nu_termwise(X); Y = replace(J, ALG); check(Y)
    assert cohom(Y) == cohom(J), "replacement is not a quasi-isomorphism"
    Y = minimise(Y); check(Y)
    return Y

def stalk(t):
    J = Cx(); J.T[0] = [t]; return J

def serre_check(X, Y):
    """dim H^t Hom(X, P_i) == dim H^{-t}(Y)_i, Y = nu X"""
    HY = cohom(Y); ks = X.degs()
    for i in range(1, NV + 1):
        Pi = ALG.proj(i)
        for t in range(-ks[-1] - 1, -ks[0] + 2):
            # C^t = Hom(X^{-t}, P_i); d: C^t -> C^{t+1} is precomposition with d^{-t-1}
            src = [a for a, s in enumerate(X.T.get(-t, [])) if std_ok(s, Pi)]
            def dmat(tt):   # C^tt -> C^{tt+1}
                s1 = [a for a, s in enumerate(X.T.get(-tt, [])) if std_ok(s, Pi)]
                s2 = [a for a, s in enumerate(X.T.get(-tt - 1, [])) if std_ok(s, Pi)]
                if not s1 or not s2: return np.zeros((len(s2), len(s1)), dtype=np.int64)
                M = X.mat(-tt - 1)   # rows deg -tt, cols deg -tt-1
                return M[np.ix_(s1, s2)].T
            h = len(src) - rank(dmat(t)) - rank(dmat(t - 1)) if src else 0
            assert h == HY.get((i, -t), 0), ("Serre mismatch", i, t, h, HY.get((i, -t), 0))

def objects(alg, spec):
    out = []
    for s in spec.split(","):
        if s == "S": out += [("S%d" % l, (l, l)) for l in range(1, alg.n + 1)]
        elif s == "P": out += [("P%d" % l, alg.proj(l)) for l in range(1, alg.n + 1)]
        elif s == "U": out += [("U%d-%d" % (a, e), (a, e)) for a in range(1, alg.n + 1) for e in range(a, alg.U[a] + 1)]
        elif s[0] == "S": l = int(s[1:]); out.append((s, (l, l)))
        elif s[0] == "P": l = int(s[1:]); out.append((s, alg.proj(l)))
        elif s[0] == "U": a, e = map(int, s[1:].split("-")); out.append((s, (a, e)))
    seen = set(); res = []
    for name, iv in out:
        if iv not in seen: seen.add(iv); res.append((name, iv))
    return res

if __name__ == "__main__":
    a, serre = args()
    rel = a[0]; ms = [int(x) for x in a[1].split(",")]; spec = a[2] if len(a) > 2 else "S,P"
    ALG = Alg(rel); NV = ALG.n
    print("REL", rel, "n", NV, "P", P, "plant", PLANT, "U", ALG.U[1:], flush=True)
    t0 = time.time(); res = {m: [] for m in ms}; sizes = {m: [] for m in ms}
    for name, iv in objects(ALG, spec):
        X = minimise(replace(stalk(iv), ALG)); check(X)
        H0 = cohom(X); want = {(v, 0): 1 for v in range(iv[0], iv[1] + 1)}
        assert H0 == want, (name, H0)
        k = 0
        for m in ms:
            while k < m:
                Y = nu(X)
                if serre: serre_check(X, Y)
                X = Y; k += 1
            H = cohom(X); ts = {t for (_, t) in H}
            sizes[m].append(sum(len(v) for v in X.T.values()))
            if len(ts) == 1 and H == {(v, list(ts)[0]): 1 for v in range(iv[0], iv[1] + 1)}:
                res[m].append((name, -list(ts)[0]))
        print(" ", name, "done", f"{time.time()-t0:.1f}s", flush=True)
    for m in ms:
        sh = sorted({c for _, c in res[m]})
        print(rel, "m", m, "objects", len(sizes[m]), "returns", len(res[m]), "shifts", sh,
              {c: [nm for nm, cc in res[m] if cc == c] for c in sh}, "max summands", max(sizes[m]), f"{time.time()-t0:.0f}s", flush=True)
