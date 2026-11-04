"""Round 115 (theorist), T10 / H-015: Hochschild cohomology dims HH^0, HH^1, HH^2 of a bound quiver algebra kQ/I (Q acyclic, finite).
HH^m is a derived invariant (Rickard/Happel); computed from the E-relative reduced bar complex Hom_{E-E}(r^{(x)m}, A) over GF(P) (P = 32003, a
generic characteristic stand-in for Q; the coefficients are the repo's presented ones).  Independent of the Coxeter key and of the Cartan matrix.
API: hh(vertices, edges, rels, top=2) -> (dim A, [dim HH^0, HH^1, ...]); edges = [(s,t,k)], rels = [{path(tuple of arrows): coef}]  (path = arrows in order)."""
import itertools, numpy as np
P = 32003
def inv(a): return pow(int(a) % P, P - 2, P)
def rank_ker(M):
    """rank of M (rows x cols) mod P, in place elimination on a copy."""
    M = np.array(M, dtype=np.int64) % P
    if M.size == 0: return 0
    r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = None
        nz = np.nonzero(M[r:, c])[0]
        if nz.size == 0: continue
        piv = r + nz[0]
        if piv != r: M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * inv(M[r, c])) % P
        idx = np.nonzero(M[:, c])[0]
        for i in idx:
            if i != r: M[i] = (M[i] - M[i, c] * M[r]) % P
        r += 1
        if r == rows: break
    return r
class Alg:
    def __init__(self, V, edges, rels):
        self.V = sorted(V); self.edges = sorted(edges)
        out = {v: [] for v in self.V}
        for e in self.edges: out[e[0]].append(e)
        paths = [((), v, v) for v in self.V]   # (arrows, s, t)
        front = [((e,), e[0], e[1]) for e in self.edges]
        while front:
            paths += front; nf = []
            for p, s, t in front:
                for e in out[t]: nf.append((p + (e,), s, e[1]))
            front = nf
            if len(paths) > 200000: raise RuntimeError('too many paths (cycle?)')
        self.paths = paths; self.pi = {p[0] if p[0] else ('e', p[1]): i for i, p in enumerate(paths)}
        def key(p, v=None): return p if p else ('e', v)
        # ideal: u r w
        gens = {}
        for rel in rels:
            first = next(iter(rel)); s = first[0][0]; t = first[-1][1]
            for pa, s0, t0 in paths:
                if pa and False: pass
        # relation endpoints
        relinfo = []
        for rel in rels:
            first = next(iter(rel)); relinfo.append((rel, first[0][0], first[-1][1]))
        byend_s = {}; byend_t = {}
        for pa, s, t in paths: byend_s.setdefault(s, []).append((pa, s, t)); byend_t.setdefault(t, []).append((pa, s, t))
        vec = []   # (block (s,t), dict path->coef)
        for rel, rs, rt in relinfo:
            for ua, us, ut in byend_t.get(rs, []):
                for wa, ws, wt in byend_s.get(rt, []):
                    d = {}
                    for pth, c in rel.items():
                        full = ua + tuple(pth) + wa; d[full] = (d.get(full, 0) + c) % P
                    vec.append(((us, wt), d))
        blocks = {}
        for pa, s, t in paths: blocks.setdefault((s, t), []).append(key(pa, s))
        self.basis = []; self.nf = {}   # nf[(s,t)] = (echelon rows, pivots, basis list)
        self.bidx = {}
        for blk, plist in blocks.items():
            col = {p: i for i, p in enumerate(plist)}
            rows = []
            for b, d in vec:
                if b == blk:
                    r = np.zeros(len(plist), dtype=np.int64)
                    for pth, c in d.items(): r[col[key(pth, blk[0])]] += c
                    rows.append(r % P)
            piv = []; E = []
            for r in rows:
                for (pc, er) in zip(piv, E):
                    if r[pc]: r = (r - r[pc] * er) % P
                nz = np.nonzero(r)[0]
                if nz.size == 0: continue
                c = nz[-1]   # pivot on the LAST nonzero so normal forms use early (short) paths
                r = (r * inv(r[c])) % P
                for k in range(len(E)):
                    if E[k][c]: E[k] = (E[k] - E[k][c] * r) % P
                piv.append(c); E.append(r)
            free = [i for i in range(len(plist)) if i not in piv]
            self.nf[blk] = (col, plist, piv, E, free)
            for i in free: self.bidx[plist[i]] = len(self.basis); self.basis.append((plist[i], blk[0], blk[1]))
        self.n = len(self.basis)
        self.rad = [i for i, (p, s, t) in enumerate(self.basis) if p[0] != 'e']
        self.keyf = key
    def reduce(self, pth, s, t):
        """vector (dict basis index -> coef) of a path in A"""
        col, plist, piv, E, free = self.nf[(s, t)]
        r = np.zeros(len(plist), dtype=np.int64); r[col[self.keyf(pth, s)]] = 1
        for pc, er in zip(piv, E):
            if r[pc]: r = (r - r[pc] * er) % P
        return {self.bidx[plist[i]]: int(r[i]) for i in free if r[i]}
    def mul(self, a, b):
        """product of basis elements a then b (concatenation), as dict; zero if endpoints mismatch."""
        pa, sa, ta = self.basis[a]; pb, sb, tb = self.basis[b]
        if ta != sb: return {}
        if pa[0] == 'e': return {b: 1}
        if pb[0] == 'e': return {a: 1}
        key = (pa + pb, sa, tb)
        if not hasattr(self, '_m'): self._m = {}
        k = (a, b)
        if k not in self._m: self._m[k] = self.reduce(pa + pb, sa, tb)
        return self._m[k]
def hh(V, edges, rels, top=2, maxdim=6000):
    A = Alg(V, edges, rels); R = A.rad
    S = lambda i: A.basis[i][1]; T = lambda i: A.basis[i][2]
    # chains of radical basis elements, composable
    chains = {0: [()]}
    chains[1] = [(r,) for r in R]
    for m in range(2, top + 2):
        chains[m] = [c + (r,) for c in chains[m - 1] for r in R if T(c[-1]) == S(r)]
    # C^m basis: (chain, target) ; chain of length 0 = a vertex e_v, target in e_v A e_v
    def cbasis(m):
        out = []
        if m == 0:
            for v in A.V:
                for b in range(A.n):
                    if S(b) == v and T(b) == v: out.append(((v,), b))
            return out
        for c in chains[m]:
            for b in range(A.n):
                if S(b) == S(c[0]) and T(b) == T(c[-1]): out.append((c, b))
        return out
    CB = {m: cbasis(m) for m in range(0, top + 2)}
    for m in CB:
        if len(CB[m]) > maxdim: return A.n, None, {m: len(CB[m]) for m in CB}
    CI = {m: {x: i for i, x in enumerate(CB[m])} for m in CB}
    def dmat(m):
        """matrix of d: C^m -> C^{m+1}: rows = basis of C^m (phi = indicator of (chain,target)), cols = C^{m+1}."""
        M = np.zeros((len(CB[m]), len(CB[m + 1])), dtype=np.int64)
        # (d phi)(a1..a_{m+1}) = a1 phi(a2..) + sum (-1)^i phi(..a_i a_{i+1}..) + (-1)^{m+1} phi(a1..am) a_{m+1}
        for cj, (chain, b0) in enumerate(CB[m + 1]):
            pass
        idx = {}
        for (c, b), i in CI[m].items(): idx.setdefault(c, []).append((b, i))
        for j, (ch, tgt) in enumerate(CB[m + 1]):
            a = ch   # a1..a_{m+1}; we need coefficient of basis cochain of C^{m+1} indexed (ch, tgt) in d(phi), for each phi basis -> loop phi
        # build by acting: for each target chain ch of length m+1 and each phi=(c,b): contributions to (ch, x)
        for ch in chains[m + 1]:
            # left term: a1 * phi(a2..a_{m+1})
            terms = []
            c = ch[1:] if m >= 1 else (S(ch[0]) if False else None)
            if m == 0:
                # phi in C^0: phi(e_v)=b; (d phi)(a1) = a1 phi(e_{t(a1)}) - phi(e_{s(a1)}) a1
                a1 = ch[0]
                for (cc, b), i in CI[0].items():
                    v = cc[0]
                    if v == T(a1):
                        for x, co in A.mul(a1, b).items(): M[i, CI[1][((a1,), x)]] += co
                    if v == S(a1):
                        for x, co in A.mul(b, a1).items(): M[i, CI[1][((a1,), x)]] -= co
                continue
            for b, i in idx.get(ch[1:], []):
                for x, co in A.mul(ch[0], b).items(): M[i, CI[m + 1][(ch, x)]] += co
            for k in range(m):   # merge a_{k+1} a_{k+2}  (0-based k,k+1), sign (-1)^{k+1}
                prod = A.mul(ch[k], ch[k + 1])
                for y, cy in prod.items():
                    if y not in set(R): continue
                    c2 = ch[:k] + (y,) + ch[k + 2:]
                    for b, i in idx.get(c2, []):
                        M[i, CI[m + 1][(ch, b)]] += ((-1) ** (k + 1)) * cy
            for b, i in idx.get(ch[:-1], []):
                for x, co in A.mul(b, ch[-1]).items(): M[i, CI[m + 1][(ch, x)]] += ((-1) ** (m + 1)) * co
        return M % P
    ranks = {}
    for m in range(0, top + 1): ranks[m] = rank_ker(dmat(m))
    dims = []
    for m in range(0, top + 1):
        ker = len(CB[m]) - ranks[m]; im = ranks[m - 1] if m else 0
        dims.append(ker - im)
    return A.n, dims, {m: len(CB[m]) for m in CB}
if __name__ == '__main__':
    # self-tests: A_n path algebra (HH = k, 0, 0); Kronecker (HH^0=1, HH^1=3); k[x]/x^2 not acyclic -> skip; commutative square
    print(hh([1, 2, 3], [(1, 2, 0), (2, 3, 0)], []))
    print(hh([1, 2], [(1, 2, 0), (1, 2, 1)], []))
    sq = [(1, 2, 0), (1, 3, 0), (2, 4, 0), (3, 4, 0)]
    print('comm square', hh([1, 2, 3, 4], sq, [{((1, 2, 0), (2, 4, 0)): 1, ((1, 3, 0), (3, 4, 0)): -1}]))
    print('square ab=0', hh([1, 2, 3, 4], sq, [{((1, 2, 0), (2, 4, 0)): 1}]))
    print('A3/rad^2', hh([1, 2, 3], [(1, 2, 0), (2, 3, 0)], [{((1, 2, 0), (2, 3, 0)): 1}]))
