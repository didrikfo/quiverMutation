"""Round 033 (toolsmith): base rate of an LNA/dual-LNA Coxeter key.
Usage: toolsmith_baserate.py layers m          -- E-123 layered family, ALL members (incl. circuit-free), key in LNA keys by category
       toolsmith_baserate.py walk n cls maxexp -- gate-admitted out-degree 2 walk rows: parent key / raw child key vs LNA keys (parent is in class by construction)
"""
import sys
sys.path.insert(0, '.'); MY = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec'))
from collections import Counter
mode = MY[1]
def lnaKeys(N):
    s = set()
    for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)): s.add(search._coxeterKeyOrNone(alg))
    return s
if mode == 'layers':
    exec(compile(open('workshop/rounds/031/theorist_layers.py').read().split("V = 1; AS")[0].split("def partitions")[1].join(["def partitions", ""]) if False else "", 'x', 'exec'))
    m = int(MY[2]); n = m + 4; LK = lnaKeys(n)
    L = open('workshop/rounds/031/theorist_layers.py').read()
    exec(compile("def partitions" + L.split("def partitions")[1].split("tab = Counter()")[0], 'tl', 'exec'))
    tab = Counter(); keys = Counter()
    for l1 in partitions(range(m)):
        for l2 in partitions(range(m)):
            sh = shape(l1, l2)
            rels = []
            for T, lab in ((T1, l1), (T2, l2)):
                blocks = {}
                for k, l in enumerate(lab): blocks.setdefault(l, []).append(k)
                for l, ks in blocks.items():
                    if l == 0: rels += [[[V, AS[k], v, T]] for k in ks]
                    else: rels += [[[V, AS[ks[j]], v, T], [V, AS[ks[j + 1]], v, T]] for j in range(len(ks) - 1)]
            A = build(arrows, rels, list(range(1, n + 1)))
            kd, mono = kerdim(A, v, procedure.relationsFrom(A))
            key = search._coxeterKeyOrNone(A)
            cat = 'free' if not sh else ('W-only' if all(s[0] == 'W' for s in sh) else 'nonW-circuit')
            tab[(cat, kd > 0, key in LK)] += 1
    print('layers m =', m, 'n =', n, 'LNA keys:', len(LK))
    for k, c in sorted(tab.items(), key=str): print(c, 'category', k[0], 'J!=0', k[1], 'key in LNA keys', k[2])
else:
    n, cls, maxexp = int(MY[2]), int(MY[3]), int(MY[4]); LK = lnaKeys(n)
    classes = {}
    for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)): classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
    order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
    seen = set(); frontier = []; tab = Counter(); nalg = 0
    def add(alg):
        k = fingerprint.canonicalKey(alg)
        if k is not None:
            if k in seen: return
            seen.add(k)
        frontier.append(alg)
    for s in classes[base]: add(s)
    while frontier and nalg < maxexp:
        cur, frontier[:] = frontier[:], []
        for alg in cur:
            if nalg >= maxexp: break
            nalg += 1
            if list(nx.simple_cycles(alg.quiver)): continue
            pk = search._coxeterKeyOrNone(alg)
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                raw = mutation.quiverMutationAtVertex(alg, v)
                bad = any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw))
                ch = None if bad else reduction.reducePathAlgebra(raw)
                ck = None if ch is None else search._coxeterKeyOrNone(ch)
                od = len(ap.arrowsOutOf(alg.quiver, v))
                tab[(od, 'parent in LNA keys', pk in LK, 'parent in class', pk == base, 'child legal', not bad, 'child in LNA keys', ck in LK if ck is not None else None, 'child in class', ck == base)] += 1
                if not bad and ck == base: add(ch)
    print('n', n, 'class', cls, 'expanded', nalg)
    for k, c in sorted(tab.items(), key=str): print(c, k)
