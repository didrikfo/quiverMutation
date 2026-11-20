"""T6 / H-017b pre-registration support (round 133, skeptic).  Does the proposed per-LNA test have power?

usage (from the repository root):
  .venv/bin/python workshop/rounds/133/skeptic_h017b.py units N         -> non-quipu LNAs at n = N, keys, duality units
  .venv/bin/python workshop/rounds/133/skeptic_h017b.py n9              -> per-unit mean defect at n = 9 from E-283's json
  .venv/bin/python workshop/rounds/133/skeptic_h017b.py power N [SEED]  -> planted-effect power curve, nulls, confound check
  .venv/bin/python workshop/rounds/133/skeptic_h017b.py collect N NAME  -> depth-5 member census of one LNA (E-283 walk;
                                                                          writes skeptic_defect_<N>_<NAME>.json here)
Unit = the duality orbit {A, A^op}: reachedQuipuAlgebras walks both orientations, so A and A^op give the same members
(n = 9: the 9 LNAs of E-283 are 4 twin pairs + 1, i.e. 5 independent units)."""
import sys, os, json, random, collections, statistics
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
MYARGV = list(sys.argv)


def name(r):
    return "".join(map(str, r))


def units(n):
    from quivermutation import coxeterTables as ct
    keys = ct.lnaKeys(n); q = ct.quipuClassLnas(n)
    out = sorted(r for r in keys if r not in q)
    cart = {r: np.array(ct.lnaCartanMatrix(n, r)) for r in out}
    J = np.eye(n, dtype = int)[::-1]
    byC = {cart[r].tobytes(): r for r in out}
    assert len(byC) == len(out)
    opp = {r: byC[(J @ cart[r].T @ J).astype(cart[r].dtype).tobytes()] for r in out}
    seen, us = set(), []
    for r in out:
        if r in seen: continue
        seen |= {r, opp[r]}
        us.append({"rep": r, "opp": opp[r], "key": keys[r], "C": cart[r]})
    return us


def features(u):
    """Maverick's registered LNA features (rank C dropped: unitriangular, always n) plus the extras used in the
    confound check (NOT regressed out)."""
    C, r = u["C"], u["rep"]
    return {"maxC": int(C.max()), "nnzC": int((C != 0).sum()),
            "nrel": sum(1 for x in r if x), "sumrel": sum(r), "minrel": min(x for x in r if x)}


def ftest(y, keys, X, perms, rng):
    """T = between-key variance of OLS residuals of y on [1, X]; p by permuting key labels among units
    (class sizes kept). Returns (T, p)."""
    A = np.column_stack([np.ones(len(y))] + ([X] if X is not None and X.size else []))
    res = y - A @ np.linalg.lstsq(A, y, rcond = None)[0]
    lab = np.array(keys)
    def T(l):
        g = [res[l == k] for k in set(l)]
        return sum(len(x) * (x.mean() - res.mean()) ** 2 for x in g) / max(1e-12, sum(((x - x.mean()) ** 2).sum() for x in g))
    t0 = T(lab); c = 0
    for _ in range(perms):
        if T(rng.permutation(lab)) >= t0 - 1e-12: c += 1
    return t0, (c + 1) / (perms + 1)


def design(us, rng, per = 3, minc = 3):
    """Pre-registered sample: classes with >= minc units; `per` units per class drawn without looking."""
    by = collections.defaultdict(list)
    for u in us: by[u["key"]].append(u)
    s = []
    for k, v in sorted(by.items(), key = lambda kv: len(kv[1])):
        if len(v) >= minc: s += [v[i] for i in rng.choice(len(v), min(per, len(v)), replace = False)]
    return s


def power(n, seed = 1):
    us = units(n); rng = np.random.default_rng(seed)
    global KID; KID = {k: i for i, k in enumerate(sorted({u["key"] for u in us}))}   # (hash(-1) == hash(-2): do not hash keys)
    REPS, PERMS = 200, 199
    for per in (3, 5):
        print("== n=%d design: %d units per key class with >= 3 units" % (n, per))
        rows = []
        for delta in (0, 0.5, 1, 1.5, 2, 3):
            rej = 0
            for _ in range(REPS):
                s = design(us, rng, per)
                ks = [KID[u["key"]] for u in s]; kl = sorted(set(ks))
                f = rng.standard_normal(len(kl)); f = (f - f.mean()) / f.std()
                fk = dict(zip(kl, f))
                X = np.array([[features(u)[c] for c in ("maxC", "nnzC")] for u in s], float)
                X = np.column_stack([X, rng.poisson(80, len(s))])   # members-per-unit covariate, unknown before the walk
                y = rng.standard_normal(len(s)) + delta * np.array([fk[k] for k in ks])
                rej += ftest(y, ks, X, PERMS, rng)[1] < 0.05
            rows.append((delta, rej / REPS)); print("  delta %.1f SD: rejection rate %.2f (N=%d, classes %d)" % (delta, rej / REPS, len(s), len(kl)))
        # null controls
        for lbl, gen in (("LNA-only random shift (truth NO)", lambda s: rng.standard_normal(len(s)) * 2),
                         ("regressed feature maxC + noise (truth NO beyond features)", lambda s: np.array([features(u)["maxC"] for u in s], float) + 0.3 * rng.standard_normal(len(s)))):
            rej = 0
            for _ in range(REPS):
                s = design(us, rng, per); ks = [KID[u["key"]] for u in s]
                X = np.array([[features(u)[c] for c in ("maxC", "nnzC")] for u in s], float)
                rej += ftest(gen(s), ks, X, PERMS, rng)[1] < 0.05
            print("  control %s: rejection rate %.2f" % (lbl, rej / REPS))
        # confound check: an LNA quantity OUTSIDE the feature list, with no planted key term
        for c in ("nrel", "sumrel", "minrel"):
            rej = 0
            for _ in range(REPS):
                s = design(us, rng, per); ks = [KID[u["key"]] for u in s]
                X = np.array([[features(u)[d] for d in ("maxC", "nnzC")] for u in s], float)
                y = np.array([features(u)[c] for u in s], float) + 0.3 * rng.standard_normal(len(s))
                rej += ftest(y, ks, X, PERMS, rng)[1] < 0.05
            print("  confound: response = source %s (+noise), key not planted: 'key effect' rate %.2f" % (c, rej / REPS))


def n9():
    rows = {}
    for fn in sorted(os.listdir(os.path.join(HERE, "..", "105"))):
        if fn.startswith("experimentalist_defect_") and fn.endswith(".json"):
            d = json.load(open(os.path.join(HERE, "..", "105", fn)))
            ds = [m["rels"] - m["cords"] for m in d["members"]]
            rows[d["lna"]] = (tuple(d["lnakey"]), len(ds), statistics.mean(ds), sorted(collections.Counter(ds).items()))
    for k, v in rows.items(): print(k, v[0][:4], "members", v[1], "mean %.3f" % v[2], v[3])
    sig = collections.defaultdict(list)
    for k, v in rows.items(): sig[(v[0], v[1], round(v[2], 9))].append(k)
    print("distinct (key, members, mean) units:", len(sig), [sorted(x) for x in sig.values()])
    big = [m for (k, c, m) in sig if k[2] == 0]; small = [m for (k, c, m) in sig if k[2] != 0]
    sd = statistics.stdev(big)
    print("class 1 unit means", ["%.3f" % x for x in big], "SD %.3f; class 2 (3033030) mean %.3f -> offset %.2f SD" % (sd, small[0], (small[0] - statistics.mean(big)) / sd))
    print("(SD of the class-1 unit means; the offset is in those units)")
    for lna, v in rows.items():   # E-303's shape effect at unit level: share of defect-1 and defect-5 members
        h = dict(v[3]); print(" shape", lna, v[0][:4], "frac d=1 %.3f  d=5 %.3f" % (h.get(1, 0) / v[1], h.get(5, 0) / v[1]))


def collect(n, nm):
    import networkx as nx
    from quivermutation import quipuRelations as qr, nakayama as nk, quipuForms as qf, invariants as inv, pathAlgebra
    def algOf(arrows, rels):
        a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, n + 1))
        for t, h in arrows: a.add_arrow(t + 1, h + 1)
        for p in rels: a.add_rel([[v + 1 for v in p]])
        return a
    alg = nk.LinearNakayamaAlgebra(n, [int(x) for x in nm])
    out = {"lna": nm, "lnakey": list(inv.coxeterKey(alg)), "members": []}
    for c, p in qr.reachedQuipuAlgebras(alg, 5).items():
        g = nx.Graph(); g.add_edges_from(c[0]); par = qf.quipuParameters(g)
        out["members"].append({"len": len(p), "cords": sum(1 for x in par[1] if x > 0), "rels": len(c[1]),
                               "key": list(inv.coxeterKey(algOf(c[0], c[1])))})
    json.dump(out, open(os.path.join(HERE, "skeptic_defect_%d_%s.json" % (n, nm)), "w"))
    ds = [m["rels"] - m["cords"] for m in out["members"]]
    print(nm, len(ds), "members, defects", sorted(collections.Counter(ds).items()))


if __name__ == "__main__":
    a = MYARGV[1:]
    if a[0] == "units":
        us = units(int(a[1])); by = collections.Counter(u["key"] for u in us)
        print("units", len(us), "self-dual", sum(u["rep"] == u["opp"] for u in us), "per key class", sorted(by.values()))
        for k, c in sorted(by.items(), key = lambda kv: kv[1]):
            print(" class size %d: %s" % (c, [name(u["rep"]) for u in us if u["key"] == k][:6]))
    elif a[0] == "n9": n9()
    elif a[0] == "power": power(int(a[1]), int(a[2]) if len(a) > 2 else 1)
    elif a[0] == "collect": collect(int(a[1]), a[2])
