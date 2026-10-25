"""Round 087 T4: H1 check with an INDEPENDENT mutation step.
usage: skeptic_homotopy.py LENGTH WIDTH [calib|run|ctl_after|ctl_seq|ctl_cox] [OP SIGN]

Mutation engine: none from the package.  Silting mutation in K^b(proj A) done by hand:
complexes of projectives over the monomial LNA (Hom(P_a,P_b) <= 1-dim, composition = masked
matrix product), Hom in the homotopy category by linear algebra mod p, minimal one-sided
add(M)-approximation, cone, Gaussian cancellation of contractible parts.  Verdict per
application: (i) every Hom(T,T[k]) = 0 for k != 0 (T tilting), (ii) End(T) has exactly the
Hom-dimension pattern of the predicted LNA, (iii) consecutive arrow composites are nonzero
exactly where the predicted paths are (so End(T) = predicted LNA as a labelled algebra).
Shared with the package: only the rule table VERIFIED_MOVES and anchorOf (the object audited).
"""
import sys, os, itertools
import numpy as np
from functools import lru_cache

P = 32003
args = sys.argv[1:]
L, W = int(args[0]), int(args[1])
MODE = args[2] if len(args) > 2 else 'run'
OP = int(args[3]) if len(args) > 3 else 0       # algebra side convention
SIGN = int(args[4]) if len(args) > 4 else 1     # +1: positive vertex = left mutation (cone of left approx)

# ---------------- LNA combinatorics (own code) ----------------
def admissible(n, r):
    rel = [(i + 1, a) for i, a in enumerate(r) if a]
    if any(a < 2 or s + a > n for s, a in rel): return False
    return all(s1 + a1 < s2 + a2 for (s1, a1), (s2, a2) in zip(rel, rel[1:]))

def enum(n):
    out = []
    def rec(i, prev_end, row):
        if i == n - 2:
            out.append(tuple(row)); return
        rec(i + 1, prev_end, row + [0])
        for a in range(2, n - i):
            if i + 1 + a > n: break
            if i + 1 + a > prev_end: rec(i + 1, i + 1 + a, row + [a])
    rec(0, 0, [])
    return out

def match(r, rule, ws):
    w, before = rule[0], rule[1]
    inside = []
    for i, a in enumerate(r):
        if not a: continue
        s = i + 1; lo, hi = s, s + a - 1
        if hi < ws or lo > ws + w - 1: continue
        if lo < ws or hi > ws + w - 1: return False
        inside.append((s - ws, a))
    return sorted(inside) == sorted(before)

def apply_rule(n, r, rule, ws):
    w, before, after, seq = rule[:4]
    new = list(r)
    for s, a in before: new[ws + s - 1] = 0
    for s, a in after:
        if new[ws + s - 1]: return None
        new[ws + s - 1] = a
    if not admissible(n, new): return None
    return tuple(new), [v + ws - 1 if v > 0 else v - ws + 1 for v in seq]

def pathNonzero(r, lo, hi):
    """path of vertices lo..hi (lo<=hi) nonzero iff no relation s..s+a inside."""
    return not any(a and s >= lo and s + a <= hi for s, a in ((i + 1, a) for i, a in enumerate(r)))

def makeAllowed(r):
    def allowed(a, b):          # P_a -> P_b nonzero
        lo, hi = (b, a) if OP == 0 else (a, b)
        return lo <= hi and pathNonzero(r, lo, hi)
    return allowed

# ---------------- linear algebra mod P ----------------
def rref(M):
    M = np.array(M, dtype=np.int64) % P
    if M.ndim == 1: M = M.reshape(1, -1)
    rows, cols = M.shape
    piv = []; r = 0
    for c in range(cols):
        if r >= rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        k = r + nz[0]
        if k != r: M[[r, k]] = M[[k, r]]
        M[r] = (M[r] * pow(int(M[r, c]), P - 2, P)) % P
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % P
        piv.append(c); r += 1
    return M[:r], piv

def nullspace(M, ncols):
    if M.shape[0] == 0: return [np.eye(ncols, dtype=np.int64)[i] for i in range(ncols)]
    R, piv = rref(M)
    free = [c for c in range(ncols) if c not in piv]
    out = []
    for f in free:
        v = np.zeros(ncols, dtype=np.int64); v[f] = 1
        for i, pc in enumerate(piv): v[pc] = (-R[i, f]) % P
        out.append(v)
    return out

def rankOf(vs):
    if len(vs) == 0: return 0
    return len(rref(np.array(vs))[1])

# ---------------- complexes ----------------
class Cx:
    _c = 0
    def __init__(self, terms, d):
        # terms: dict k -> list of vertices ; d: dict k -> matrix (len(terms[k+1]) x len(terms[k]))
        self.terms = {k: list(v) for k, v in terms.items() if len(v)}
        self.d = {}
        for k in self.terms:
            if k + 1 in self.terms:
                m = np.array(d.get(k, np.zeros((len(self.terms[k + 1]), len(self.terms[k])))), dtype=np.int64) % P
                self.d[k] = m.reshape(len(self.terms[k + 1]), len(self.terms[k]))
        Cx._c += 1; self.id = Cx._c
    def degs(self): return sorted(self.terms)
    def T(self, k): return self.terms.get(k, [])
    def D(self, k):
        if k in self.d: return self.d[k]
        return np.zeros((len(self.T(k + 1)), len(self.T(k))), dtype=np.int64)

def maskOf(allowed, src, tgt):
    return np.array([[1 if allowed(a, b) else 0 for a in src] for b in tgt], dtype=np.int64).reshape(len(tgt), len(src))

def mmul(allowed, G, srcF, midlist, tgtG, F):
    """G (tgt x mid) after F (mid x src), masked."""
    if G.shape[1] == 0 or F.shape[1] == 0 or G.shape[0] == 0:
        return np.zeros((G.shape[0], F.shape[1]), dtype=np.int64)
    return ((G @ F) % P) * maskOf(allowed, srcF, tgtG) % P

class World:
    def __init__(self, r):
        self.r = r; self.allowed = makeAllowed(r); self.cache = {}
    # degree-k maps X -> Y: dict m -> matrix len(Y.T(m+k)) x len(X.T(m))
    def vars(self, X, Y, k):
        out = []
        for m in X.degs():
            tg = Y.T(m + k)
            for ci, a in enumerate(X.T(m)):
                for ri, b in enumerate(tg):
                    if self.allowed(a, b): out.append((m, ri, ci))
        return out
    def vecToMap(self, X, Y, k, vs, vec):
        mp = {m: np.zeros((len(Y.T(m + k)), len(X.T(m))), dtype=np.int64) for m in X.degs() if len(Y.T(m + k))}
        for (m, ri, ci), x in zip(vs, vec): mp[m][ri, ci] = x
        return mp
    def mapToVec(self, X, Y, k, vs, mp):
        return np.array([int(mp[m][ri, ci]) if m in mp else 0 for (m, ri, ci) in vs], dtype=np.int64)
    def diff(self, X, Y, k, mp):
        """D f = d_Y f - (-1)^k f d_X , degree k+1."""
        res = {}
        for m in set(X.degs()) | set(Y.degs()):
            # component X^m -> Y^{m+k+1}
            tg = Y.T(m + k + 1); sc = X.T(m)
            if not tg or not sc: continue
            acc = np.zeros((len(tg), len(sc)), dtype=np.int64)
            if (m in mp) and Y.T(m + k) and (m + k) in Y.d:
                acc += mmul(self.allowed, Y.d[m + k], sc, Y.T(m + k), tg, mp[m])
            if (m + 1) in mp and m in X.d and X.T(m + 1):
                acc -= (1 if k % 2 == 0 else -1) * mmul(self.allowed, mp[m + 1], sc, X.T(m + 1), tg, X.d[m])
            res[m] = acc % P
        return res
    def hom(self, X, Y, k=0):
        key = (X.id, Y.id, k)
        if key in self.cache: return self.cache[key]
        vs = self.vars(X, Y, k)
        vs1 = self.vars(X, Y, k + 1)
        idx1 = {v: i for i, v in enumerate(vs1)}
        # matrix of D : maps^k -> maps^{k+1}
        cols = []
        for i, v in enumerate(vs):
            e = np.zeros(len(vs), dtype=np.int64); e[i] = 1
            mp = self.vecToMap(X, Y, k, vs, e)
            dm = self.diff(X, Y, k, mp)
            cols.append(self.mapToVec(X, Y, k + 1, vs1, dm) if vs1 else np.zeros(0, dtype=np.int64))
        Dm = np.array(cols, dtype=np.int64).T.reshape(len(vs1), len(vs)) if vs else np.zeros((len(vs1), 0), dtype=np.int64)
        Z = nullspace(Dm, len(vs)) if len(vs) else []
        # coboundaries: D of degree k-1 maps
        vsm = self.vars(X, Y, k - 1); Bs = []
        for i, v in enumerate(vsm):
            e = np.zeros(len(vsm), dtype=np.int64); e[i] = 1
            dm = self.diff(X, Y, k - 1, self.vecToMap(X, Y, k - 1, vsm, e))
            Bs.append(self.mapToVec(X, Y, k, vs, dm))
        rb = rankOf(Bs) if (Bs and len(vs)) else 0
        dimH = len(Z) - rb
        res = dict(vs=vs, Z=Z, B=Bs, dim=dimH, X=X, Y=Y, k=k)
        self.cache[key] = res
        return res
    def reps(self, X, Y):
        """chain-map representatives (dicts) of a basis of Hom_K(X,Y) complementing coboundaries."""
        h = self.hom(X, Y, 0)
        if h['dim'] == 0: return []
        cur = [b for b in h['B']] if h['B'] else []
        rk = rankOf(cur) if cur else 0
        out = []
        for z in h['Z']:
            t = cur + [z]
            r2 = rankOf(t)
            if r2 > rk: cur = t; rk = r2; out.append(self.vecToMap(X, Y, 0, h['vs'], z))
        return out
    def compose(self, X, Y, Zc, f, g):
        """g: Y->Z after f: X->Y (degree 0 maps)."""
        res = {}
        for m in X.degs():
            if m in f and m in g and Zc.T(m) and Y.T(m):
                res[m] = mmul(self.allowed, g[m], X.T(m), Y.T(m), Zc.T(m), f[m])
        return res
    def inB(self, X, Y, f):
        h = self.hom(X, Y, 0)
        v = self.mapToVec(X, Y, 0, h['vs'], f)
        if not v.any(): return True
        if not h['B']: return False
        return rankOf(h['B'] + [v]) == rankOf(h['B'])

def minimalize(allowed, terms, d):
    terms = {k: list(v) for k, v in terms.items()}
    d = {k: np.array(m, dtype=np.int64) % P for k, m in d.items()}
    changed = True
    while changed:
        changed = False
        for k in sorted(terms):
            if k not in d or k + 1 not in terms: continue
            M = d[k]
            for ri, b in enumerate(terms[k + 1]):
                for ci, a in enumerate(terms[k]):
                    if a == b and M[ri, ci] % P:
                        c0inv = pow(int(M[ri, ci]), P - 2, P)
                        rr = [i for i in range(M.shape[0]) if i != ri]; cc = [j for j in range(M.shape[1]) if j != ci]
                        gam = M[np.ix_(rr, [ci])]; dl = M[np.ix_([ri], cc)]
                        eps = M[np.ix_(rr, cc)]
                        corr = (gam @ dl) % P * c0inv % P
                        srcv = [terms[k][j] for j in cc]; tgv = [terms[k + 1][i] for i in rr]
                        corr = corr * maskOf(allowed, srcv, tgv) % P
                        d[k] = (eps - corr) % P
                        if k - 1 in d: d[k - 1] = np.delete(d[k - 1], ci, axis=0)
                        if k + 1 in d: d[k + 1] = np.delete(d[k + 1], ri, axis=1)
                        terms[k] = srcv; terms[k + 1] = tgv
                        changed = True
                        break
                if changed: break
            if changed: break
    return {k: v for k, v in terms.items() if v}, d

def checkD2(allowed, terms, d):
    for k in terms:
        if k in d and k + 1 in d and terms.get(k + 2):
            prod = mmul(allowed, d[k + 1], terms[k], terms[k + 1], terms[k + 2], d[k])
            if prod.any(): return False
    return True

class Illegal(Exception): pass

def mutate(W_, T, v, left):
    """T: dict vertex -> Cx.  Returns new Cx for position v (left: cone of left approx; right: cone(g)[-1])."""
    X = T[v]
    others = [j for j in T if j != v]
    gens = []  # (j, chain map)
    if left:
        hv = {j: W_.hom(X, T[j], 0)['dim'] for j in others}
        for j in others:
            if hv[j] == 0: continue
            if W_.hom(T[j], T[j], 0)['dim'] != 1: raise Illegal('nonlocal End')
            h = W_.hom(X, T[j], 0)
            S = list(h['B'])
            for k in others:
                if k == j or hv[k] == 0 or W_.hom(T[k], T[j], 0)['dim'] == 0: continue
                for f in W_.reps(X, T[k]):
                    for g in W_.reps(T[k], T[j]):
                        S.append(W_.mapToVec(X, T[j], 0, h['vs'], W_.compose(X, T[k], T[j], f, g)))
            rk = rankOf(S) if S else 0
            for z in h['Z']:
                t = S + [z]; r2 = rankOf(t)
                if r2 > rk: S = t; rk = r2; gens.append((j, W_.vecToMap(X, T[j], 0, h['vs'], z)))
    else:
        hv = {j: W_.hom(T[j], X, 0)['dim'] for j in others}
        for j in others:
            if hv[j] == 0: continue
            if W_.hom(T[j], T[j], 0)['dim'] != 1: raise Illegal('nonlocal End')
            h = W_.hom(T[j], X, 0)
            S = list(h['B'])
            for k in others:
                if k == j or hv[k] == 0 or W_.hom(T[j], T[k], 0)['dim'] == 0: continue
                for f in W_.reps(T[j], T[k]):
                    for g in W_.reps(T[k], X):
                        S.append(W_.mapToVec(T[j], X, 0, h['vs'], W_.compose(T[j], T[k], X, f, g)))
            rk = rankOf(S) if S else 0
            for z in h['Z']:
                t = S + [z]; r2 = rankOf(t)
                if r2 > rk: S = t; rk = r2; gens.append((j, W_.vecToMap(T[j], X, 0, h['vs'], z)))
    # E = sum of copies
    degsE = set()
    for j, _ in gens: degsE |= set(T[j].degs())
    Eterms = {m: [] for m in degsE}; Eoff = []  # offsets per gen per degree
    for j, _ in gens:
        off = {}
        for m in T[j].degs():
            off[m] = len(Eterms.setdefault(m, [])); Eterms[m] += T[j].T(m)
        Eoff.append(off)
    Ed = {}
    for m in Eterms:
        if m + 1 in Eterms:
            M = np.zeros((len(Eterms[m + 1]), len(Eterms[m])), dtype=np.int64)
            for (j, _), off in zip(gens, Eoff):
                if m in off and (m + 1) in off and m in T[j].d:
                    dj = T[j].d[m]; M[off[m + 1]:off[m + 1] + dj.shape[0], off[m]:off[m] + dj.shape[1]] = dj
            Ed[m] = M
    terms = {}; d = {}
    if left:
        # cone(f)^k = X^{k+1} + E^k ; d(t,e) = (-d_X t, f t + d_E e)
        ks = set(m - 1 for m in X.degs()) | set(Eterms)
        for k in ks: terms[k] = list(X.T(k + 1)) + list(Eterms.get(k, []))
        for k in ks:
            if k + 1 not in terms: continue
            nx0, nx1 = len(X.T(k + 1)), len(X.T(k + 2))
            ne0, ne1 = len(Eterms.get(k, [])), len(Eterms.get(k + 1, []))
            M = np.zeros((nx1 + ne1, nx0 + ne0), dtype=np.int64)
            if (k + 1) in X.d: M[:nx1, :nx0] = -X.d[k + 1]
            if k in Ed and ne0 and ne1: M[nx1:, nx0:] = Ed[k]
            # f: X^{k+1} -> E^{k+1}
            for (j, mp), off in zip(gens, Eoff):
                if (k + 1) in mp and (k + 1) in off:
                    blk = mp[k + 1]; M[nx1 + off[k + 1]: nx1 + off[k + 1] + blk.shape[0], :nx0] += blk
            d[k] = M % P
    else:
        # cone(g)[-1]^k = E^k + X^{k-1}; d = [[-d_E, 0],[g, d_X]]
        ks = set(Eterms) | set(m + 1 for m in X.degs())
        for k in ks: terms[k] = list(Eterms.get(k, [])) + list(X.T(k - 1))
        for k in ks:
            if k + 1 not in terms: continue
            ne0, ne1 = len(Eterms.get(k, [])), len(Eterms.get(k + 1, []))
            nx0, nx1 = len(X.T(k - 1)), len(X.T(k))
            M = np.zeros((ne1 + nx1, ne0 + nx0), dtype=np.int64)
            if k in Ed and ne0 and ne1: M[:ne1, :ne0] = -Ed[k]
            if (k - 1) in X.d and nx0 and nx1: M[ne1:, ne0:] = X.d[k - 1]
            for (j, mp), off in zip(gens, Eoff):
                if k in mp and k in off:
                    blk = mp[k]   # E^k copy -> X^k ; target rows in X-part of degree k+1 (= X^k)
                    M[ne1:, off[k]:off[k] + blk.shape[1]] += blk
            d[k] = M % P
    allowed = W_.allowed
    if not checkD2(allowed, terms, d): raise Illegal('d^2 != 0')
    terms2, d2 = minimalize(allowed, terms, d)
    if not checkD2(allowed, terms2, d2): raise Illegal('d^2 != 0 after cancel')
    return Cx(terms2, d2)

def initial(n):
    return {i: Cx({0: [i]}, {}) for i in range(1, n + 1)}

SPAN = [0]
def verdict(n, r, T, pred, W_, info=None):
    """None if T is tilting and End(T) = LNA(pred) as an algebra (vertices matched along the arrow line)."""
    for a in T:
        for b in T:
            lo = min(T[b].terms) - max(T[a].terms); hi = max(T[b].terms) - min(T[a].terms)   # Hom(X,Y[k])=0 outside [lo,hi]
            SPAN[0] = max(SPAN[0], hi, -lo)
            for k in range(min(-4, lo), max(4, hi) + 1):
                if k == 0: continue
                if W_.hom(T[a], T[b], k)['dim']: return 'NOTTILTING Hom(T%d,T%d[%d])!=0' % (a, b, k)
    dim = {(a, b): W_.hom(T[a], T[b], 0)['dim'] for a in T for b in T}
    if any(dim[(a, a)] != 1 for a in T): return 'NONLOCAL'
    if any(dim[(a, b)] > 1 for a in T for b in T): return 'DIM>1'
    # irreducible morphisms: Hom(a,b) not generated by composites through other summands
    arrows = {}
    for a in T:
        for b in T:
            if a == b or not dim[(a, b)]: continue
            irr = True
            for c in T:
                if c in (a, b) or not dim[(a, c)] or not dim[(c, b)]: continue
                for f in W_.reps(T[a], T[c]):
                    for g in W_.reps(T[c], T[b]):
                        if not W_.inB(T[a], T[b], W_.compose(T[a], T[c], T[b], f, g)): irr = False
            if irr: arrows[(a, b)] = W_.reps(T[a], T[b])[0]
    outd = {a: [b for (x, b) in arrows if x == a] for a in T}; ind = {a: [x for (x, b) in arrows if b == a] for a in T}
    if any(len(v) > 1 for v in outd.values()) or any(len(v) > 1 for v in ind.values()) or len(arrows) != n - 1:
        return 'NOTLINE quiver arrows=%s' % sorted(arrows)
    starts = [a for a in T if not ind[a]]
    if len(starts) != 1: return 'NOTLINE'
    line = [starts[0]]
    while outd[line[-1]]: line.append(outd[line[-1]][0])
    if len(line) != n: return 'NOTLINE'
    pos = {}
    for idx, a in enumerate(line): pos[(n - idx) if OP == 0 else (idx + 1)] = a   # position -> vertex
    if info is not None: info['perm'] = tuple(pos[i] for i in range(1, n + 1))
    allowedP = makeAllowed(pred)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if dim[(pos[i], pos[j])] != (1 if allowedP(i, j) else 0): return 'DIM (%d,%d) got %d want %d' % (i, j, dim[(pos[i], pos[j])], 1 if allowedP(i, j) else 0)
    # composites of consecutive arrows along the line
    for s0 in range(n - 1):
        chain = arrows[(line[s0], line[s0 + 1])]
        for e in range(s0 + 2, n):
            chain = W_.compose(T[line[s0]], T[line[e - 1]], T[line[e]], chain, arrows[(line[e - 1], line[e])])
            nz = not W_.inB(T[line[s0]], T[line[e]], chain)
            pi, pj = [k for k, v in pos.items() if v == line[s0]][0], [k for k, v in pos.items() if v == line[e]][0]
            if nz != allowedP(pi, pj): return 'COMPOSITE %d->%d nonzero=%s want=%s' % (pi, pj, nz, allowedP(pi, pj))
    return None

def run_one(n, r, seq, pred, left_for_positive, info=None):
    W_ = World(r)
    T = initial(n)
    for s in seq:
        left = (s > 0) == left_for_positive
        T[abs(s)] = mutate(W_, T, abs(s), left)
    return verdict(n, r, T, pred, W_, info)

if __name__ == '__main__':
    from quivermutation import lnaMoves as lm
    rows = enum(L)
    rules = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == W]
    leftPos = (SIGN == 1)
    tot = 0; cnt = {}; first = None; examples = {}
    import time; t0 = time.time()
    # Coxeter-equal alternatives for ctl_cox (own Cartan, exact via sympy)
    if MODE == 'ctl_cox':
        import sympy as sp
        @lru_cache(None)
        def cox(r):
            n = L; rel = [(i + 1, a) for i, a in enumerate(r) if a]
            C = sp.zeros(n, n)
            for i in range(1, n + 1):
                for j in range(i, n + 1):
                    if not any(i <= s and s + a <= j for s, a in rel): C[i - 1, j - 1] = 1
            Phi = -C.T * C.inv()
            return tuple(sp.Poly(Phi.charpoly(sp.Symbol('x')).as_expr(), sp.Symbol('x')).all_coeffs())
        byCox = {}
        for q in rows: byCox.setdefault(cox(q), []).append(q)
    for ri, rule in enumerate(rules):
        for r in rows:
            for ws in range(1, L - W + 1):
                if not match(r, rule, ws): continue
                ap = apply_rule(L, r, rule, ws)
                if ap is None: continue
                pred, seq = ap
                if MODE == 'ctl_after':
                    nxt = rules[(ri + 1) % len(rules)]
                    pred2 = apply_rule(L, r, (rule[0], rule[1], nxt[2], rule[3]), ws)
                    if pred2 is None or pred2[0] == pred: continue
                    pred = pred2[0]
                elif MODE == 'ctl_seq':
                    seq = [(s + 1 if s > 0 else s - 1) if abs(s) < L else s for s in seq]
                elif MODE == 'ctl_cox':
                    alts = [q for q in byCox[cox(pred)] if q != pred]
                    if not alts: continue
                    pred = alts[(hash(pred) % 7919) % len(alts)]
                try:
                    v = run_one(L, r, seq, pred, leftPos)
                except Illegal as e:
                    v = 'ILLEGAL ' + str(e)
                tot += 1
                key = 'OK' if v is None else v.split(' ')[0]
                cnt[key] = cnt.get(key, 0) + 1
                if v is not None and first is None: first = (rule, r, ws, v)
                if v is not None and cnt[key] <= 3: examples.setdefault(key, []).append((rule[1], rule[2], rule[3], r, ws, v))
    print('MODE', MODE, 'L', L, 'W', W, 'OP', OP, 'SIGN', SIGN, 'rules', len(rules), 'apps', tot, cnt, 'time %.0fs' % (time.time() - t0), 'max|k| needed', SPAN[0])
    for k, v in examples.items():
        for e in v[:3]: print('  ex', k, e)
