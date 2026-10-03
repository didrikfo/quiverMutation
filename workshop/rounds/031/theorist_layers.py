"""Round 031 (theorist): exhaustive family  i -> {a_1..a_m} -> v -> {t1, t2}  (n = m+4) with, for each sink t_j, an arbitrary set partition of the m paths i a_k v t_j
plus a 'zero' block (monomial + commutativity relations, scalar 1).  Gamma_i is read off the partitions; J = kerdim from the library; key vs LNA keys.
Question: does any member with a circuit other than the length-2 ground path (W) -- an nn 2-cycle or a longer circuit -- have a Coxeter key in an LNA/dual-LNA key set of length n?
Usage: theorist_layers.py m"""
import sys, itertools
sys.path.insert(0, '.'); MY = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec'))
from collections import Counter
m = int(MY[1]); npend = int(MY[2]) if len(MY) > 2 else 0; n = m + 4
lk = set()
lkN = {}
for N in range(n, n + npend + 1):
    lkN[N] = set()
    for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)): lkN[N].add(search._coxeterKeyOrNone(alg))
def partitions(items):  # set partitions of items plus optional zero block: yield label list with label 0 = zero
    items = list(items)
    def rec(i, lab, used):
        if i == len(items): yield tuple(lab); return
        for l in range(0, used + 2):
            if l == used + 1: yield from rec(i + 1, lab + [l], used + 1)
            else: yield from rec(i + 1, lab + [l], used)
    yield from rec(0, [], 0)
V = 1; AS = list(range(2, 2 + m)); v = 2 + m; T1, T2 = v + 1, v + 2
arrows = [(V, a) for a in AS] + [(a, v) for a in AS] + [(v, T1), (v, T2)]
def shape(l1, l2):
    # vertices ('1',l) l>0, ('2',l) l>0, ground 'g'; edge k joins node(l1[k]) and node(l2[k])
    nd = lambda s, l: 'g' if l == 0 else (s, l)
    par = {}
    def f(x):
        par.setdefault(x, x)
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    edges = [(nd('1', a), nd('2', b)) for a, b in zip(l1, l2)]
    for x, y in edges: par[f(x)] = f(y)
    comp = {}
    for x, y in edges: comp.setdefault(f(x), []).append((x, y))
    res = []
    for c, es in comp.items():
        vs = {z for e in es for z in e}
        if len(es) >= len(vs):   # has a circuit (cyclomatic number >= 1)
            kind = 'W' if (len(es) == 2 and 'g' in vs and len(vs) == 2) else ('nn' if len(es) == 2 else 'long%d' % len(es))
            res.append((kind, len(es), len(vs)))
    return tuple(sorted(res))
tab = Counter(); hits = []; cnt = 0
for l1 in partitions(range(m)):
    for l2 in partitions(range(m)):
        sh = shape(l1, l2)
        if not sh: continue
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
        inl = key in lkN[n]
        if kd > 0 and not mono and npend:
            # pendant extensions (not at v): any LNA key at the longer length?
            def ext(arr_, nodes_, depth):
                if depth == npend: return False
                new = max(nodes_) + 1; got = False
                for u in nodes_:
                    if u == v: continue
                    for e in ((u, new), (new, u)):
                        B = build(arr_ + [e], rels, nodes_ + [new])
                        if search._coxeterKeyOrNone(B) in lkN[len(nodes_) + 1]: got = True
                        got = ext(arr_ + [e], nodes_ + [new], depth + 1) or got
                return got
            inl = ext(arrows, list(range(1, n + 1)), 0) or inl
        tab[(sh, kd > 0, mono, inl)] += 1
        if inl and not (sh and all(s[0] == 'W' for s in sh)) and not mono: hits.append((l1, l2, sh))
print('m =', m, 'n =', n)
for k, c in sorted(tab.items(), key=str): print(c, 'shape', k[0], 'J!=0', k[1], 'single-path-in-J', k[2], 'key in LNA keys', k[3])
print('non-W circuit with LNA key and no single path in J:', hits[:10], len(hits))
