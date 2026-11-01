"""H-017 / T6: is the defect (relations - cords) constant per Coxeter-key class, n = 9, depth 5?
usage: collect NAME   -> forward walk of LNA NAME to depth 5, writes experimentalist_defect_<NAME>.json
       analyse        -> census + controls from the saved json files (and E-203/E-221 member txt files)
Members = monomial quipu algebras reached (quipuRelations.reachedQuipuAlgebras), any path length, cords >= 0."""
import sys, json, glob, os, random, collections, itertools
import networkx as nx
HERE = os.path.dirname(os.path.abspath(__file__))
N, DEPTH = 9, 5
NAMES = ["3033030","3345000","3455000","3505000","4444400","4550400","5040330","5050030","5504030"]

def algOf(arrows, rels):
    from quivermutation import pathAlgebra
    a = pathAlgebra.PathAlgebra(); a.add_vertices_from(range(1, N + 1))
    for t, h in arrows: a.add_arrow(t + 1, h + 1)
    for p in rels: a.add_rel([[v + 1 for v in p]])
    return a

def collect(name):
    from quivermutation import quipuRelations as qr, nakayama as nk, quipuForms as qf, invariants as inv
    alg = nk.LinearNakayamaAlgebra(N, [int(x) for x in name])
    lnakey = list(inv.coxeterKey(alg))
    out = {"lna": name, "lnakey": lnakey, "members": []}
    for c, p in qr.reachedQuipuAlgebras(alg, DEPTH).items():
        g = nx.Graph(); g.add_edges_from(c[0]); par = qf.quipuParameters(g)
        cords = sum(1 for x in par[1] if x > 0); rels = len(c[1])
        key = list(inv.coxeterKey(algOf(c[0], c[1])))
        out["members"].append({"len": len(p), "cords": cords, "rels": rels, "key": key})
    json.dump(out, open(os.path.join(HERE, "experimentalist_defect_%s.json" % name), "w"))
    print(name, len(out["members"]), "members")

def load():
    return {n: json.load(open(os.path.join(HERE, "experimentalist_defect_%s.json" % n))) for n in NAMES
            if os.path.exists(os.path.join(HERE, "experimentalist_defect_%s.json" % n))}

def stat(rows):
    """number of key classes (key-tuple) in which defect is not constant; rows = [(key, defect)]"""
    d = collections.defaultdict(set)
    for k, x in rows: d[k].add(x)
    return sum(1 for s in d.values() if len(s) > 1), len(d)

def analyse():
    data = load(); rows = []
    for n, d in data.items():
        for m in d["members"]:
            if m["key"] != d["lnakey"]: print("KEY MISMATCH", n, m)
            rows.append((tuple(m["key"]), m["rels"] - m["cords"], n, m["cords"], m["len"]))
    print("members", len(rows))
    cls = collections.defaultdict(list)
    for r in rows: cls[r[0]].append(r)
    print("key classes", len(cls), "LNAs per class", {k[:4]: sorted({x[2] for x in v}) for k, v in cls.items()})
    for k, v in cls.items():
        print(" class", k[:4], "n=", len(v), "defects", sorted(collections.Counter(x[1] for x in v).items()))
    for n in data:
        v = [r for r in rows if r[2] == n]
        print(" LNA", n, "defect counts", sorted(collections.Counter(x[1] for x in v).items()),
              "min", min(x[1] for x in v), "min with cords>=1", min([x[1] for x in v if x[3] >= 1] or [None]))
    obs = stat([(r[0], r[1]) for r in rows])
    print("OBS: classes with non-constant defect / classes:", obs)
    # filtered view mimicking E-203/E-221: cords>=1, len==5, defect<=1
    flt = [(r[0], r[1]) for r in rows if r[3] >= 1 and r[4] == 5 and r[1] <= 1]
    print("E-203/221-style filter (cords>=1, L=5, defect<=1):", len(flt), "members, defects", set(x[1] for x in flt), "non-const classes", stat(flt))
    random.seed(1)
    # control D (planted error): in the constant filtered set, flip ONE member's defect to 2; the test must now say NO
    f2 = list(flt); f2[0] = (f2[0][0], 2)
    print("planted single defect error in filtered set:", stat(f2), "(expect 1 non-constant)")
    # tallies per LNA of the filtered set, for comparison with E-203/E-221 member counts
    print("filtered members per LNA:", sorted(collections.Counter(r[2] for r in rows if r[3] >= 1 and r[4] == 5 and r[1] <= 1).items()))
    # control A: shuffled key labels. Planted truth: defect is NOT a function of key unless it is.
    # test: does the statistic distinguish real labels from shuffled ones?  Also a planted function of key.
    base = [(r[0], r[1]) for r in rows]
    keys = [k for k, _ in base]; ds = [x for _, x in base]
    sh = []
    for _ in range(1000):
        random.shuffle(keys); sh.append(stat(list(zip(keys, ds)))[0])
    print("shuffled labels: non-constant classes mean %.2f min %d max %d (obs %d)" % (sum(sh)/len(sh), min(sh), max(sh), obs[0]))
    # control B: planted constancy -- replace defect by a function of key; test must say constant (YES)
    ids = {k: i for i, k in enumerate(sorted(set(keys)))}
    pl = [(r[0], ids[r[0]] % 3) for r in rows]
    print("planted defect = f(key):", stat(pl), "(expect 0 non-constant)")
    # control C: planted class split -- take the key classes and split each by a hidden label (cords parity of member);
    # statistic on the SPLIT classes must change the verdict if the hidden label carries the defect
    spl = [((r[0], r[3]), r[1]) for r in rows]
    print("classes split by cords:", stat(spl))
    spl2 = [((r[0], r[1]), r[1]) for r in rows]
    print("classes split by defect itself (trivially constant):", stat(spl2))

if __name__ == "__main__":
    if sys.argv[1] == "collect": collect(sys.argv[2])
    else: analyse()
