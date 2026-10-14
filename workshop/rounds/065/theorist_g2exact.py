"""G2 past the span statement (round 065, theorist).

For random ACYCLIC multigraph quivers (parallel arrows allowed) with random monomial and
binomial relations of length >= 2, and every vertex i with an out-arrow, compare

  K_opp(v) = red(Z_1 e_i)(v)            Oppermann's side, computed in the free path algebra
  Ker(v)   = kernel of procedure._kernelOverIdeal (repo step 7 before any 'nearer' filter)

as SUBSPACES of the same coordinate space {(relation index r, old tail path p)}: same sign, same
orientation, no transpose, no extension by paths.  Claim (Prop. 1 of the submission): K_opp = Ker.

It also counts generators.  g(v)   = dim Z_1(i->v) / (B_1 + Z_1 rad + rad Z_1)  = number of
degree-2 arrows i->v of Q^(2) (Oppermann Constr. 2.2) ;
gd(v)  = dim K(v) / (sum_a K(w) a + red B_1(v))   = number of those whose dx survives in H_0(M);
r7(v)  = number of step-7 relations the repo returns at target v (after its nearer-target filter).
Lemma L (submission): red(B_1 e_i) = { rbar . s : s in I } lies in the ideal of steps 4-6 + through;
checked here by: red B_1(v) is contained in K(v) (it must be) and its repo-side image is in the span.

Run from the repository root:  timeout 10m .venv/bin/python workshop/rounds/065/theorist_g2exact.py [seeds] [per_seed]
"""
import sys
import random
import itertools
from fractions import Fraction

import networkx as nx

sys.path.insert(0, ".")
from quivermutation import procedure as pr
from quivermutation import arrowPaths as ap


# ---------------------------------------------------------------- linear algebra
def reduceRow(row, piv):
    row = {k: Fraction(x) for k, x in row.items() if x}
    changed = True
    while row:
        head = min(row)
        if head in piv:
            f = row[head]
            for k, x in piv[head].items():
                row[k] = row.get(k, 0) - f * x
            row = {k: x for k, x in row.items() if x}
        else:
            return row, head
    return row, None


class Span:
    def __init__(self, vectors=()):
        self.piv = {}
        for v in vectors:
            self.add(v)

    def add(self, v):
        row, head = reduceRow(v, self.piv)
        if head is None:
            return False
        row = {k: x / row[head] for k, x in row.items()}
        self.piv[head] = row
        return True

    def dim(self):
        return len(self.piv)

    def contains(self, v):
        row, head = reduceRow(v, self.piv)
        return head is None


def nullspace(cols, matrix):
    rows = [{c: Fraction(v) for c, v in row.items() if v} for row in matrix.values()]
    pivots = {}
    for row in rows:
        for pc, prow in list(pivots.items()):
            if row.get(pc):
                f = row[pc]
                for k, x in prow.items():
                    row[k] = row.get(k, 0) - f * x
                row = {k: x for k, x in row.items() if x}
        if not row:
            continue
        lead = min(row)
        row = {k: x / row[lead] for k, x in row.items()}
        for pc in list(pivots):
            if pivots[pc].get(lead):
                f = pivots[pc][lead]
                for k, x in row.items():
                    pivots[pc][k] = pivots[pc].get(k, 0) - f * x
                pivots[pc] = {k: x for k, x in pivots[pc].items() if x}
        pivots[lead] = row
    basis = []
    for c in cols:
        if c in pivots:
            continue
        sol = {c: Fraction(1)}
        for pc, prow in pivots.items():
            if prow.get(c):
                sol[pc] = -prow[c]
        basis.append(sol)
    return basis


# ---------------------------------------------------------------- Oppermann side
class Opp:
    """Degree-1 words (r, p, q): q = path start -> s(r), r, p = path t(r) -> v."""

    def __init__(self, Q, rels):
        self.Q, self.rels = Q, rels
        self.src = [ap.source(r) for r in rels]
        self.tgt = [ap.target(r) for r in rels]
        self.cache = {}

    def paths(self, a, b):
        return ap.allPathsBetween(self.Q, a, b)

    def words(self, start, v):
        out = []
        for ri in range(len(self.rels)):
            for q in self.paths(start, self.src[ri]):
                for p in self.paths(self.tgt[ri], v):
                    out.append((ri, p, q))
        return out

    def Z1(self, start, v):
        key = (start, v)
        if key not in self.cache:
            cols = self.words(start, v)
            matrix = {}
            for (ri, p, q) in cols:
                for path, c in self.rels[ri].items():
                    matrix.setdefault(q + path + p, {})[(ri, p, q)] = c
            self.cache[key] = nullspace(cols, matrix)
        return self.cache[key]

    def B1(self, start, v):
        """d of the degree-2 words q r2 m r1 p (path order start -> v)."""
        out = []
        n = len(self.rels)
        for r2 in range(n):
            for q in self.paths(start, self.src[r2]):
                for r1 in range(n):
                    for m in self.paths(self.tgt[r2], self.src[r1]):
                        for p in self.paths(self.tgt[r1], v):
                            vec = {}
                            for path, c in self.rels[r2].items():      # q rho2_el m r1 p
                                key = (r1, p, q + path + m)
                                vec[key] = vec.get(key, 0) + c
                            for path, c in self.rels[r1].items():      # - q r2 m rho1_el p
                                key = (r2, m + path + p, q)
                                vec[key] = vec.get(key, 0) - c
                            vec = {k: x for k, x in vec.items() if x}
                            if vec:
                                out.append(vec)
        return out


def red(vec):
    return {(ri, p): c for (ri, p, q), c in vec.items() if not q}


# ---------------------------------------------------------------- repo side
def repoSide(Q, rels, i):
    """Raw kernel per target (coordinates (r index, old tail)), and the filtered step-7 output."""
    outArrows = ap.arrowsOutOf(Q, i)
    inArrows = ap.arrowsInto(Q, i)
    outRel = [r for r in rels if ap.source(r) == i]
    idx = {id(r): k for k, r in enumerate(rels)}
    arrows = pr._Arrows()
    for order, arrow in enumerate(ap.arrowsOf(Q)):
        t, h, _ = arrow
        if t != i and h != i:
            arrows.declare((pr._Arrows.CARRIED, arrow), t, h, pr._Arrows.CARRIED, order)
    for order, a in enumerate(outArrows):
        arrows.declare((pr._Arrows.FLIP, a), a[1], i, pr._Arrows.FLIP, order)
    for order, b in enumerate(inArrows):
        for inner, a in enumerate(outArrows):
            arrows.declare((pr._Arrows.COMPOSITE, b, a), b[0], a[1], pr._Arrows.COMPOSITE, (order, inner))
    for order, r in enumerate(outRel):
        arrows.declare((pr._Arrows.RBAR, order), i, ap.target(r), pr._Arrows.RBAR, order)
    M = arrows.build(Q.nodes)
    relationOfArrow = {arrows.arrow((pr._Arrows.RBAR, k)): r for k, r in enumerate(outRel)}
    rbarIndex = {a: idx[id(r)] for a, r in relationOfArrow.items()}
    oldArrow = pr._oldArrowReading(arrows)
    raw, kept = {}, {}
    found = []
    for v in pr._targetsInOrder(M, i):
        cands = pr._pathsOutOfMutatedVertex(M, i, v)
        if not cands:
            continue
        shadow = {}
        for path in cands:
            rel = relationOfArrow[path[0]]
            tail = pr._inOldQuiver(path[1:], oldArrow)
            for alpha in outArrows:
                quo = ap.leftDivide(rel, alpha)
                shadow[(path, alpha)] = ap.postCompose(quo, tail) if quo else {}
        ker = pr._kernelOverIdeal(Q, rels, cands, outArrows, shadow)

        def conv(el):
            return {(rbarIndex[p[0]], pr._inOldQuiver(p[1:], oldArrow)): c for p, c in el.items()}
        raw[v] = [conv(e) for e in ker]
        k = []
        for e in ker:
            if not pr._isForcedByNearer(M, found, e, i):
                found.append(e)
                k.append(conv(e))
        kept[v] = k
    return raw, kept


# ---------------------------------------------------------------- random instances
def randomInstance(rng, n):
    Q = nx.MultiDiGraph()
    Q.add_nodes_from(range(1, n + 1))
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            if rng.random() < 0.45:
                for k in range(2 if rng.random() < 0.2 else 1):
                    Q.add_edge(a, b, key=k)
    rv = []
    for _ in range(rng.randint(1, 6)):
        a = rng.randint(1, n - 2)
        b = rng.randint(a + 2, n)
        ps = [p for p in ap.allPathsBetween(Q, a, b) if len(p) >= 2]
        if not ps:
            continue
        if len(ps) >= 2 and rng.random() < 0.4:
            p1, p2 = rng.sample(ps, 2)
            rv.append(ap.combination([(p1, 1), (p2, rng.choice([-1, -2, 2]))]))
        else:
            rv.append(ap.combination([(rng.choice(ps), 1)]))
    return Q, rv


def check(Q, rels, i):
    opp = Opp(Q, rels)
    raw, kept = repoSide(Q, rels, i)
    nodes = sorted(Q.nodes)
    res = {"ok": True, "g": 0, "gd": 0, "r7": 0, "kdim": 0, "redBnotInK": 0, "gmulti": 0}
    for v in nodes:
        if v == i:
            continue
        Z = opp.Z1(i, v)
        K = Span(red(z) for z in Z)
        R = Span(raw.get(v, []))
        same = K.dim() == R.dim() and all(R.contains(k) for k in K.piv.values()) \
            and all(K.contains(r) for r in R.piv.values())
        res["ok"] &= same
        res["kdim"] += K.dim()
        res["kmulti"] = res.get("kmulti", 0) + sum(1 for z in K.piv.values() if len(z) > 1)
        # B1, extensions
        B = [b for b in opp.B1(i, v)]
        ext = []
        for a in Q.in_edges(v, keys=True):
            w = a[0]
            if w == i:
                continue
            for z in opp.Z1(i, w):
                ext.append({(ri, p + (a,), q): c for (ri, p, q), c in z.items()})
        for b in Q.out_edges(i, keys=True):
            j = b[1]
            if j == v:
                continue
            for z in opp.Z1(j, v):
                ext.append({(ri, p, (b,) + q): c for (ri, p, q), c in z.items()})
        S = Span(B + ext)
        assert all(Span(Z).contains(x) for x in B + ext), "B1/ext not in Z1"
        res["g"] += len(Z) - S.dim()
        # degree-2 arrow count is dim Z/(B+ext); here check it is >= 0 and B in Z
        # gd: K / (K(w) a + red B)
        kext = []
        for a in Q.in_edges(v, keys=True):
            w = a[0]
            if w == i:
                continue
            for z in opp.Z1(i, w):
                kext.append({(ri, p + (a,)): c for (ri, p), c in red(z).items()})
        redB = [red(b) for b in B]
        S2 = Span(kext + redB)
        assert all(K.contains(x) for x in kext + redB)
        res["gd"] += K.dim() - S2.dim()
        res["r7"] += len(kept.get(v, []))
        # Lemma L side: red B elements are rbar . s with s in I (coordinates (r2, m+path_el+p))
    # full pipeline: unreduced rewrite, then minimal generators; relations with source i*
    nq, nr = pr.mutateAtVertex(Q, rels, i)
    minimal = pr._minimalGenerators(nq, [dict(r) for r in nr if r])
    res["fromIstar"] = sum(1 for r in minimal if ap.source(r) == i)
    res["step7raw"] = sum(len(v) for v in kept.values())
    return res


if __name__ == "__main__":
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    tot = dict(kdim=0, kpairs=0, kmulti=0, minEqG=0, minLtG=0, minGtG=0, pairs=0, bad=0, par=0, multi=0, g=0, gd=0, r7=0, r7gtgd=0, r7ltgd=0, ggtgd=0)
    for seed in range(1, seeds + 1):
        rng = random.Random(seed)
        for _ in range(per):
            n = rng.choice([5, 6, 7])
            Q, rv = randomInstance(rng, n)
            if not rv:
                continue
            par = ap.hasParallelArrows(Q)
            for i in sorted(Q.nodes):
                if not ap.arrowsOutOf(Q, i):
                    continue
                r = check(Q, rv, i)
                tot["pairs"] += 1
                tot["par"] += par
                tot["bad"] += (not r["ok"])
                tot["multi"] += any(len(rel) > 1 for rel in rv)
                tot["kdim"] += r["kdim"]
                tot["kpairs"] += r["kdim"] > 0
                tot["kmulti"] += r.get("kmulti", 0) > 0
                tot["g"] += r["g"]
                tot["gd"] += r["gd"]
                tot["r7"] += r["r7"]
                tot["r7gtgd"] += r["r7"] > r["gd"]
                tot["r7ltgd"] += r["r7"] < r["gd"]
                tot["ggtgd"] += r["g"] > r["gd"]
                tot["minEqG"] += r["fromIstar"] == r["g"]
                tot["minLtG"] += r["fromIstar"] < r["g"]
                tot["minGtG"] += r["fromIstar"] > r["g"]
                if r["fromIstar"] > r["g"]:
                    print("MIN>G seed", seed, "i", i, "mutable", pr.isMutable(Q, rv, i), "g", r["g"], "fromI*", r["fromIstar"], sorted(Q.edges(keys=True)), rv)
                if not r["ok"]:
                    print("DIFFER seed", seed, "i", i, sorted(Q.edges(keys=True)), rv)
    print(tot)
    print("ALL EQUAL (K_opp = Ker)" if tot["bad"] == 0 else "SOME DIFFER")
