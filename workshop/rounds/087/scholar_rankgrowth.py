"""Round 087 scholar: rank growth of the minimal projective complexes nu^k P_l for a linear Nakayama algebra (E-231 proposal).
nu = Serre functor = - (x)^L D(Lambda).  A complex of projectives X is sent termwise to a complex of injectives, which is
resolved by the (functorial) E-relative normalised bar resolution and minimised by Gaussian elimination over GF(p).
Total rank R_k = sum_l (number of indecomposable summands of nu^k P_l).  Exact periodicity nu^a P ~ P[b] shows as bounded R_k.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/087/scholar_rankgrowth.py RELSTRING [K] [--cap N]
   RELSTRING = digits of the relation lengths (class name tail, e.g. 90000000 for n = 10, 0 = no relation); n = len + 2."""
import sys, collections
P = 2147483647
def kupisch(rl):
    n = len(rl) + 2
    e = []
    for i in range(1, n + 1):
        m = n
        for j in range(i, n + 1):
            if j - 1 < len(rl) and rl[j - 1] > 0: m = min(m, j + rl[j - 1] - 1)
        e.append(m)
    return n, [e[i] - (i + 1) + 1 for i in range(n)]       # p[i-1] = length of P_i
class Alg:
    def __init__(self, rl):
        self.n, self.p = kupisch(rl); n = self.n
        self.p = [None] + self.p                           # 1-indexed
        self.s = [None] + [min(s for s in range(1, n + 1) if s + self.p[s] > a) for a in range(1, n + 1)]  # I_a = [s_a, a]
        self._ch = {}
    def hom(self, a, b): return b <= a < b + self.p[b]
    def chains(self, last):
        """all tuples (x0 > ... > xs = last), consecutive Hom(x_{i-1}, x_i) != 0"""
        if last in self._ch: return self._ch[last]
        out = [(last,)]
        for y in range(last + 1, min(self.n, last + self.p[last] - 1) + 1):   # y -> last needs last <= y < last+p[last]
            for c in self.chains_up(y, last): pass
        return None
def build_chains(A):
    # chains with fixed LAST vertex, extended on the left (x0 largest)
    n = A.n; ch = {}
    for last in range(n, 0, -1):
        out = [(last,)]
        for y in range(last + 1, min(n, last + A.p[last] - 1) + 1):
            for c in ch[y]:                                  # c starts at y... chain (x0..y) then append last: need hom(y,last)
                pass
        ch[last] = out
    return ch
def chains_ending(A, last, memo):
    if last in memo: return memo[last]
    # chains (x0 > ... > x_{s-1} > last): choose x_{s-1} = y > last with hom(y,last); prefix = chain ending at y
    out = [(last,)]
    for y in range(last + 1, A.n + 1):
        if A.hom(y, last):
            for c in chains_ending(A, y, memo): out.append(c + (last,))
    memo[last] = out
    return out
class Cx:
    """complex of projectives: verts[u] = (deg, vertex); out[u] = {v: coeff}; inn[v] = {u: coeff}"""
    def __init__(self): self.vert = {}; self.out = collections.defaultdict(dict); self.inn = collections.defaultdict(dict)
    def add(self, u, deg, vtx): self.vert[u] = (deg, vtx)
    def setd(self, u, v, c):
        c %= P
        if c: self.out[u][v] = (self.out[u].get(v, 0) + c) % P; self.inn[v][u] = self.out[u][v]
        if self.out[u].get(v, 1) == 0: del self.out[u][v]; del self.inn[v][u]
def minimise(A, X):
    viol = 0
    work = [u for u in X.vert]
    alive = set(X.vert)
    stack = list(X.vert)
    while stack:
        u = stack.pop()
        if u not in alive: continue
        found = None
        for v, c in X.out[u].items():
            if X.vert[v][1] == X.vert[u][1]: found = (v, c); break
        if not found: continue
        v, c = found; ci = pow(c, P - 2, P)
        ws = [(w, cw) for w, cw in X.inn[v].items() if w != u]
        zs = [(z, cz) for z, cz in X.out[u].items() if z != v]
        for w, cw in ws:
            for z, cz in zs:
                if not A.hom(X.vert[w][1], X.vert[z][1]):
                    viol += 1; continue
                nv = (X.out[w].get(z, 0) - cw * cz % P * ci) % P
                if nv: X.out[w][z] = nv; X.inn[z][w] = nv
                else:
                    X.out[w].pop(z, None); X.inn[z].pop(w, None)
            stack.append(w)
        for x in (u, v):
            for z in list(X.out[x]): X.inn[z].pop(x, None)
            for w in list(X.inn[x]): X.out[w].pop(x, None)
            X.out.pop(x, None); X.inn.pop(x, None); alive.discard(x)
        for z, _ in zs: stack.append(z)
    Y = Cx()
    for u in alive: Y.vert[u] = X.vert[u]
    for u in alive:
        for v, c in X.out.get(u, {}).items(): Y.out[u][v] = c; Y.inn[v][u] = c
    return Y, viol
def nu_step(A, X, memo, cap):
    """apply nu then resolve + minimise"""
    T = Cx(); idx = {}
    def sid(*key):
        if key not in idx: idx[key] = len(idx)
        return idx[key]
    for u, (t, a) in X.vert.items():
        for xs in range(A.s[a], a + 1):
            for ch in chains_ending(A, xs, memo):
                T.add(sid(u, ch), t - (len(ch) - 1), ch[0])
    if len(T.vert) > cap: return None, len(T.vert), 0
    for (u, ch), k in list(idx.items()):
        t, a = X.vert[u]; s = len(ch) - 1
        if s >= 1: T.setd(k, idx[(u, ch[1:])], 1)                     # d_0 : P_{x0} -> P_{x1}
        for i in range(1, s):                                           # d_i : drop x_i if composite nonzero
            if A.hom(ch[i - 1], ch[i + 1]): T.setd(k, idx[(u, ch[:i] + ch[i + 1:])], (-1) ** i)
        if s >= 1 and ch[-2] <= a: T.setd(k, idx[(u, ch[:-1])], (-1) ** s)             # d_s : drop last (module action)
        for v, c in X.out.get(u, {}).items():                           # nu(h): I_a -> I_b, g(x) = 1 iff x <= b
            if ch[-1] <= X.vert[v][1]: T.setd(k, idx[(v, ch)], c * (-1) ** s)
    Y, viol = minimise(A, T)
    ren = {}
    Z = Cx()
    for u, dv in Y.vert.items(): ren[u] = len(ren); Z.vert[ren[u]] = dv
    for u in Y.vert:
        for v, c in Y.out[u].items(): Z.out[ren[u]][ren[v]] = c; Z.inn[ren[v]][ren[u]] = c
    return Z, len(T.vert), viol
def d2check(X):
    for u in X.vert:
        acc = collections.defaultdict(int)
        for v, c in X.out.get(u, {}).items():
            for w, c2 in X.out.get(v, {}).items(): acc[w] += c * c2
        # only meaningful where Hom allowed; caller filters
    return True
def run(rl, K, cap=400000, verbose=True):
    A = Alg(rl); memo = {}; res = []
    cur = {}
    for l in range(1, A.n + 1):
        X = Cx(); X.vert[0] = (0, l); cur[l] = X
    res.append([len(X.vert) for X in cur.values()])
    for k in range(1, K + 1):
        nxt = {}; big = 0; viol = 0
        for l, X in cur.items():
            Z, size, v = nu_step(A, X, memo, cap)
            if Z is None: nxt = None; break
            nxt[l] = Z; big = max(big, size); viol += v
        if nxt is None: res.append(None); break
        cur = nxt
        ranks = [len(X.vert) for X in cur.values()]
        spans = [max(d for d, _ in X.vert.values()) - min(d for d, _ in X.vert.values()) + 1 for X in cur.values()]
        res.append(ranks)
        if verbose: print(f"k={k}: total rank {sum(ranks)}  max {max(ranks)}  max deg-span {max(spans)}  biggest resolved complex {big}  masked violations {viol}", flush=True)
    return A, res, cur
if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    rl = [int(c) for c in args[0]]; K = int(args[1]) if len(args) > 1 else 8
    cap = int(sys.argv[sys.argv.index("--cap") + 1]) if "--cap" in sys.argv else 400000
    A, res, cur = run(rl, K, cap)
    print("n", A.n, "Kupisch p", A.p[1:])
    print("rank totals", [sum(r) if r else None for r in res])
