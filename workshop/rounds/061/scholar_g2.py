"""G2 by direct computation: repo step 7 versus the d of Oppermann's degree-2 generators.

For a quiver Q with relations I (acyclic, no loops) and a vertex i, Oppermann
(arXiv:1504.02617, Sect. 2 Construction 2.3 + Sect. 3) takes the dg quiver with a
degree-1 arrow r per minimal relation (d r = rho_r) and degree-2 arrows x that
generate the degree-1 cycles Z_1.  In the mutated quiver M a degree-2 arrow
x: i -> v has ||x|| = 1 and  partial x = dec(red(d x)) = sum_P eps_P (s_P rbar_P)
(Delta never occurs: acyclic, a path out of i has no arrow entering i first).

Oppermann side, computed here from scratch:  Z_1 e_i  =  ker of
   d(p r q) = p rho_r q   on degree-1 words  p r q  (i --q--> s(r), r, t(r) --p--> v)
in the FREE path algebra, and  K_opp(v) = red(Z_1 e_i) = {sum eps s r with q empty}.
Repo side: procedure._relationsOutOfMutatedVertex, read back into old paths.
Claim tested:  span of repo step-7 relations (extended on the left by paths) = K_opp(v),
coefficient by coefficient (same sign, same orientation, no transpose).

Run from the repository root:  .venv/bin/python workshop/rounds/061/scholar_g2.py
"""
import sys
import itertools
from fractions import Fraction

import networkx as nx

sys.path.insert(0, ".")
from quivermutation import procedure as pr
from quivermutation import arrowPaths as ap


def build(edges, relations_vertex):
    """edges: list of (t, h); relations: list of dict {vertex-tuple: coeff}."""
    Q = nx.MultiDiGraph()
    for t, h in edges:
        Q.add_edge(t, h)
    rels = []
    for rel in relations_vertex:
        comb = {}
        for vs, c in rel.items():
            comb[tuple((vs[k], vs[k + 1], 0) for k in range(len(vs) - 1))] = c
        rels.append(ap.combination(comb))
    return Q, rels


def rank(vectors):
    rows = [dict(v) for v in vectors]
    piv = {}
    r = 0
    for row in rows:
        row = {k: Fraction(x) for k, x in row.items() if x}
        while row:
            head = min(row)
            if head in piv:
                f = row[head]
                for k, x in piv[head].items():
                    row[k] = row.get(k, 0) - f * x
                row = {k: x for k, x in row.items() if x}
            else:
                row = {k: x / row[head] for k, x in row.items()}
                piv[head] = row
                r += 1
                break
    return r


def nullspace(cols, matrix):
    """matrix: dict row -> dict col -> val. Returns basis of kernel as dicts col->val."""
    rows = [dict((c, Fraction(v)) for c, v in row.items() if v) for row in matrix.values()]
    pivots = {}
    for row in rows:
        for pc, prow in pivots.items():
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


def oppermannK(Q, rels, i):
    """{v: list of vectors keyed (r_index, tail_path)}  =  red(Z_1 e_i)."""
    out = {}
    for v in sorted(Q.nodes):
        if v == i:
            continue
        cols = []
        matrix = {}
        for ri, rho in enumerate(rels):
            s, t = ap.source(rho), ap.target(rho)
            for q in ap.allPathsBetween(Q, i, s):
                for p in ap.allPathsBetween(Q, t, v):
                    cols.append((ri, p, q))
                    for path, c in rho.items():
                        matrix.setdefault(q + path + p, {})[(ri, p, q)] = c
        ker = nullspace(cols, matrix)
        red = []
        for z in ker:
            vec = {(ri, p): c for (ri, p, q), c in z.items() if q == ()}
            if vec:
                red.append(vec)
        out[v] = red
    return out


def repoStep7(Q, rels, i):
    """Step 7 of the repo, read back as {v: vectors keyed (r_index, old tail)}."""
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
    found = pr._relationsOutOfMutatedVertex(Q, rels, i, M, outArrows, outRel, arrows)
    rbarIndex = {arrows.arrow((pr._Arrows.RBAR, k)): idx[id(r)] for k, r in enumerate(outRel)}
    reading = pr._oldArrowReading(arrows)
    res = []
    for rel in found:
        vec = {}
        for path, c in rel.items():
            vec[(rbarIndex[path[0]], pr._inOldQuiver(path[1:], reading))] = c
        res.append((ap.target(rel), vec))
    return res


def compare(name, Q, rels, i):
    opp = oppermannK(Q, rels, i)
    repo = repoStep7(Q, rels, i)
    # repo relations extended on the left (post-composed) by old paths to every v
    ext = {v: [] for v in opp}
    for w, vec in repo:
        for v in opp:
            for u in ap.allPathsBetween(Q, w, v):
                ext[v].append({(ri, s + u): c for (ri, s), c in vec.items()})
    ok = True
    detail = []
    for v in opp:
        a, b = rank(opp[v]), rank(ext[v])
        u = rank(opp[v] + ext[v])
        same = (a == b == u)
        ok &= same
        if a or b:
            detail.append((v, a, b, u))
    nrepo = len(repo)
    print("%-34s i=%d  repo step-7 gens=%d  per-target (v, dimK_opp, dimK_repo, dim union): %s  %s"
          % (name, i, nrepo, detail, "EQUAL" if ok else "DIFFER"))
    return ok, repo


if __name__ == "__main__":
    allok = True
    cases = []
    lin = [(1, 2), (2, 3), (3, 4), (4, 5)]
    cases.append(("A5 rel 123, 234", lin, [{(1, 2, 3): 1}, {(2, 3, 4): 1}]))
    cases.append(("A5 rel 123, 234, 345", lin, [{(1, 2, 3): 1}, {(2, 3, 4): 1}, {(3, 4, 5): 1}]))
    cases.append(("A5 rel 1234, 2345 (len 3)", lin, [{(1, 2, 3, 4): 1}, {(2, 3, 4, 5): 1}]))
    cases.append(("A5 rel 123, 2345", lin, [{(1, 2, 3): 1}, {(2, 3, 4, 5): 1}]))
    cases.append(("D5 (Oppermann Ex 4.2)", [(1, 2), (2, 3), (3, 4), (2, 5)], [{(1, 2, 3): 1}, {(2, 3, 4): 1}]))
    cases.append(("square+tail, commutativity", [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)],
                  [{(1, 2, 4): 1, (1, 3, 4): -1}, {(2, 4, 5): 1}]))
    cases.append(("square+tail, comm + 345", [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)],
                  [{(1, 2, 4): 1, (1, 3, 4): -1}, {(2, 4, 5): 1}, {(3, 4, 5): 1}]))
    for name, edges, rv in cases:
        Q, rels = build(edges, rv)
        for i in sorted(Q.nodes):
            ok, _ = compare(name, Q, rels, i)
            allok &= ok
    import random
    bad = multi = total = 0
    for seed in (1, 2):
        rng = random.Random(seed)
        for _ in range(100):
            n = rng.choice([5, 6])
            edges = [(a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1) if rng.random() < 0.45]
            Q0 = nx.MultiDiGraph(); Q0.add_nodes_from(range(1, n + 1)); Q0.add_edges_from(edges)
            rv = []
            for _r in range(rng.randint(1, 4)):
                a = rng.randint(1, n - 2); b = rng.randint(a + 2, n)
                ps = [p for p in nx.all_simple_paths(Q0, a, b) if len(p) >= 3]
                if not ps:
                    continue
                if len(ps) >= 2 and rng.random() < 0.4:
                    p1, p2 = rng.sample(ps, 2); rv.append({tuple(p1): 1, tuple(p2): rng.choice([-1, -2, 2])})
                else:
                    rv.append({tuple(rng.choice(ps)): 1})
            if not rv:
                continue
            Q, rels = build(edges, rv)
            for i in sorted(Q.nodes):
                ok, repo = compare("random seed %d" % seed, Q, rels, i) if False else (None, None)
                opp = oppermannK(Q, rels, i); repo = repoStep7(Q, rels, i)
                ext = {v: [] for v in opp}
                for w, vec in repo:
                    for v in opp:
                        for u in ap.allPathsBetween(Q, w, v):
                            ext[v].append({(ri, s + u): c for (ri, s), c in vec.items()})
                total += 1
                multi += any(len(vec) > 1 for _w, vec in repo)
                if any(not (rank(opp[v]) == rank(ext[v]) == rank(opp[v] + ext[v])) for v in opp):
                    bad += 1
    print("random: %d pairs, %d with a multi-term step-7 generator, %d DIFFER" % (total, multi, bad))
    allok &= bad == 0
    print("ALL EQUAL" if allok else "SOME DIFFER")
