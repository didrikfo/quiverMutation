"""Round 131, theorist: thinness lemma.  If the gate admits v and
d_i = dim e_i A e_v <= 1 for every i != v, then every J_i = 0 (tilting).

Independent engine (exact Fractions, no library maths) on random acyclic
kQ/I with zero and commutativity relations, compared with the library's
gate and isTiltingAtVertex; planted failures for each dropped hypothesis;
cyclic monomial check of the path-simplicity step.

usage: theorist_thin.py sweep [N] [seed] | cyclic [N] [seed] | planted
"""
import random, sys, warnings
from fractions import Fraction
from itertools import product
import networkx as nx

# ---------- independent engine: arrows (tail, head, key), paths = tuples of arrows
def paths(arrows, i, j, maxlen=12):
    out, stack = [], [((), i)]
    while stack:
        p, x = stack.pop()
        if x == j and p: out.append(p)
        if len(p) >= maxlen: continue
        for a in arrows:
            if a[0] == x: stack.append((p + (a,), a[1]))
    return out

def src(p): return p[0][0]
def tgt(p): return p[-1][1]

class Alg:
    """Acyclic kQ/I, I = two-sided ideal of the given combinations (dict path->coef)."""
    def __init__(self, n, arrows, rels):
        self.n, self.arrows, self.rels = n, arrows, rels
        self._piv = {}
    def P(self, i, j):
        return paths(self.arrows, i, j)
    def pivots(self, i, j):
        if (i, j) in self._piv: return self._piv[(i, j)]
        span = []
        for r in self.rels:
            a, b = src(next(iter(r))), tgt(next(iter(r)))
            for u in ([()] + self.P(i, a) if i != a else [()]) if True else []:
                if u and src(u) != i: continue
                if not u and a != i: continue
                for w in ([()] + self.P(b, j)) if b != j else [()]:
                    if w and tgt(w) != j: continue
                    if not w and b != j: continue
                    span.append({u + p + w: c for p, c in r.items()})
        piv = {}
        for v in span:
            v = self._red(v, piv)
            if v:
                h = min(v); c = v[h]
                v = {k: x / c for k, x in v.items()}
                for k2 in list(piv):           # keep fully reduced
                    if h in piv[k2]:
                        f = piv[k2][h]
                        piv[k2] = {k: piv[k2].get(k, 0) - f * v.get(k, 0) for k in set(piv[k2]) | set(v)}
                        piv[k2] = {k: x for k, x in piv[k2].items() if x}
                piv[h] = v
        self._piv[(i, j)] = piv
        return piv
    @staticmethod
    def _red(v, piv):
        v = {k: Fraction(x) for k, x in v.items() if x}
        for h in sorted(piv):
            if h in v:
                f = v[h]
                v = {k: v.get(k, 0) - f * piv[h].get(k, 0) for k in set(v) | set(piv[h])}
                v = {k: x for k, x in v.items() if x}
        return v
    def nf(self, comb):
        p = next(iter(comb)); return self._red(comb, self.pivots(src(p), tgt(p)))
    def std(self, i, j):
        piv = self.pivots(i, j); return [p for p in self.P(i, j) if p not in piv]
    def out(self, v): return [a for a in self.arrows if a[0] == v]
    def gate(self, v):
        if not self.out(v): return False
        for i in range(self.n):
            if i == v: continue
            for p in self.P(i, v):
                if self.nf({p: 1}) and not any(self.nf({p + (b,): 1}) for b in self.out(v)):
                    return False
        return True
    def dJ(self, i, v):
        S = self.std(i, v)
        rows = []
        for s in S:
            row = {}
            for b in self.out(v):
                for k, x in self.nf({s + (b,): 1}).items(): row[(b, k)] = x
            rows.append(row)
        piv, rank = {}, 0
        for r in rows:
            r = self._red(r, piv)
            if r:
                h = min(r); piv[h] = {k: x / r[h] for k, x in r.items()}; rank += 1
        return len(S), len(S) - rank

# ---------- library face
def toLib(alg):
    from quivermutation import pathAlgebra, arrowPaths as ap
    pa = pathAlgebra.PathAlgebra()
    pa.quiver = nx.MultiDiGraph(); pa.quiver.add_nodes_from(range(1, alg.n + 1))
    keymap = {}
    for a in alg.arrows:
        keymap[a] = (a[0] + 1, a[1] + 1, pa.quiver.add_edge(a[0] + 1, a[1] + 1))
    pa.arrowRels = [{tuple(keymap[a] for a in p): c for p, c in r.items()} for r in alg.rels]
    pa.rels = [[list(ap.projectPath(p)) for p in r] for r in pa.arrowRels]
    return pa

def libGateTilt(pa, v):
    from quivermutation import mutation
    return mutation.mutationIsPossibleAtVertex(pa, v + 1), mutation.isTiltingAtVertex(pa, v + 1)

# ---------- random acyclic algebras
def randomAlg(rng, monomial=False):
    n = rng.randint(4, 7)
    arrows, cnt = [], {}
    for _ in range(rng.randint(n - 1, n + 3)):
        a = rng.randrange(n - 1); b = rng.randint(a + 1, min(n - 1, a + 3))
        k = cnt.get((a, b), 0)
        if k >= 2: continue
        cnt[(a, b)] = k + 1; arrows.append((a, b, k))
    alg = Alg(n, arrows, [])
    allp = [p for i in range(n) for j in range(n) if i < j for p in alg.P(i, j) if len(p) >= 2]
    rels = []
    for _ in range(rng.randint(0, 3)):
        if not allp: break
        p = rng.choice(allp)
        par = [q for q in allp if q != p and src(q) == src(p) and tgt(q) == tgt(p)]
        if par and not monomial and rng.random() < 0.6:
            rels.append({p: 1, rng.choice(par): rng.choice([-1, 1, 2])})
        else:
            rels.append({p: 1})
    return Alg(n, arrows, rels)

def sweep(N, seed, monomial=False):
    rng = random.Random(seed)
    table, mism, l1viol = {}, 0, 0
    for _ in range(N):
        alg = randomAlg(rng, monomial)
        pa = toLib(alg)
        for v in range(alg.n):
            if not alg.out(v): continue
            g = alg.gate(v)
            ds = [alg.dJ(i, v) for i in range(alg.n) if i != v]
            dmax = max(d for d, _ in ds); J = sum(j for _, j in ds)
            if g and any(j > max(d - 1, 0) for d, j in ds): l1viol += 1
            lg, lt = libGateTilt(pa, v)
            if lg != g or lt != (J == 0): mism += 1
            key = (min(dmax, 2), g, J > 0)
            table[key] = table.get(key, 0) + 1
    return table, mism, l1viol

# ---------- cyclic monomial: path-simplicity step
def cyclic(N, seed):
    """Monomial finite-dim kQ/I with cycles: nonzero paths = no relation subword.
    Check: if d_w <= 1 for all w != v (all nonzero paths counted), every nonzero
    path into v from i != v is vertex-simple.  Also count the converse failure."""
    rng = random.Random(seed)
    tested = held = nonsimple_when_thick = 0
    for _ in range(N):
        n = rng.randint(3, 5)
        arrows = []
        for _ in range(rng.randint(n, n + 2)):
            a, b = rng.randrange(n), rng.randrange(n)
            if a != b: arrows.append((a, b, len(arrows)))
        # zero relations: all paths of length 4 + a few of length 2,3 -> finite dim
        def zero(p, rel):
            return any(p[k:k + len(r)] == r for r in rel for k in range(len(p) - len(r) + 1))
        short = [p for i in range(n) for j in range(n) for p in paths(arrows, i, j, 3) if len(p) in (2, 3)]
        rel = set(rng.sample(short, min(len(short), rng.randint(0, 4))))
        L = 4
        nz = lambda i, j: [p for p in paths(arrows, i, j, L - 1) if not zero(p, rel)]
        for v in range(n):
            d = {w: len(nz(w, v)) for w in range(n) if w != v}
            ps = [p for w in d for p in nz(w, v)]
            simple = all(len({a[0] for a in p} | {p[-1][1]}) == len(p) + 1 for p in ps)
            if max(d.values(), default=0) <= 1:
                tested += 1; held += simple
            elif not simple:
                nonsimple_when_thick += 1
    return tested, held, nonsimple_when_thick

# ---------- planted
def planted():
    out = []
    # (a) E-080 square: d = 2, gate admits, J = 1
    A = [(0, 1, 0), (0, 2, 0), (1, 3, 0), (2, 3, 0), (3, 4, 0)]
    alg = Alg(5, A, [{(A[0], A[2], A[4]): 1, (A[1], A[3], A[4]): -1}])
    out.append(("E-080 square v=3 (d<=1 dropped)", alg, 3))
    # (b) gate dropped: 0->1->2, ab = 0, v = 1: d = 1, gate refuses, J = 1
    B = [(0, 1, 0), (1, 2, 0)]
    out.append(("1->2->3, ab=0, v=2 (gate dropped)", Alg(3, B, [{(B[0], B[1]): 1}]), 1))
    # (c) monomial, d = 2: gate <=> J = 0 still
    alg = Alg(5, A, [{(A[0], A[2], A[4]): 1}])
    out.append(("square, monomial abd e = 0, v=3", alg, 3))
    for name, alg, v in out:
        ds = [alg.dJ(i, v) for i in range(alg.n) if i != v]
        print(f"{name}: gate={alg.gate(v)} dmax={max(d for d,_ in ds)} J={sum(j for _,j in ds)} lib={libGateTilt(toLib(alg), v)}")
    # (d) E-318 parent, library only
    sys.path.insert(0, 'tests')
    import test_search_soundness as t
    from quivermutation import mutation, arrowPaths as ap, procedure
    P = t._e318Parent(); q = P.quiver; rel = procedure.relationsFrom(P)
    for v in sorted(q.nodes):
        if not mutation.mutationIsPossibleAtVertex(P, v): continue
        d = max([len(ap.allPathsBetween(q, i, v)) - len(ap.idealBasis(q, rel, i, v))
                 for i in q.nodes if i != v and ap.allPathsBetween(q, i, v)] or [0])
        print(f"E-318 parent v={v}: gate=True d={d} tilting={mutation.isTiltingAtVertex(P, v)}")

if __name__ == "__main__":
    warnings.simplefilter("ignore")
    mode = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 131
    if mode == "sweep":
        for mono in (False, True):
            t, m, l1 = sweep(N, seed, mono)
            print(f"monomial={mono}: (min(dmax,2), gate, J!=0) -> count {dict(sorted(t.items()))}; "
                  f"library mismatches {m}; E-128 L1 violations {l1}")
    elif mode == "cyclic":
        print("cyclic monomial: (thin vertices tested, path-simplicity held, thick vertices with a nonzero non-simple path) =",
              cyclic(N, seed))
    else:
        planted()
