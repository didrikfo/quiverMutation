"""Class labels for every LNA of length n (n <= 10): Coxeter key + known-equivalence orbits (free, edges, doubles).
A label is the key, except that a key holding several orbits is split by orbit only when the orbits contain a quipu seed of different quipus
(cospectral case); otherwise reported as ambiguous. usage: python workshop/rounds/027/maverick_classes.py N"""
import sys, collections
from quivermutation import coxeterTables as ct, freeMoves
def labels(n):
    keys = ct.lnaKeys(n)
    lnas, orbits = freeMoves.derivedOrbits(n, rules=[], free=True, edges=True, doubles=True)
    bykey = collections.defaultdict(list)
    for o, mem in orbits.items():
        bykey[keys[mem[0]]].append(mem)
        assert len({keys[m] for m in mem}) == 1
    lab = {}; amb = []
    for k, orbs in bykey.items():
        if len(orbs) > 1: amb.append((k, [len(o) for o in orbs]))
        for i, o in enumerate(orbs):
            for m in o: lab[m] = (k, i)
    return lab, amb, keys

def classes(n):
    """class label per LNA: the Coxeter key, except that a cospectral key is split into its quipu classes by the orbit
    (free+edges+doubles) joined with its mirror; returns (label dict, list of unresolved cospectral orbits)."""
    from quivermutation import quipuForms as qf, nakayama as nk
    keys = ct.lnaKeys(n)
    lnas, orbits = freeMoves.derivedOrbits(n, rules=None, free=True, edges=True, doubles=True)
    oid = {m: o for o, mem in orbits.items() for m in mem}
    par = {o: o for o in orbits}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for m in lnas:
        mr = freeMoves.mirrorRow(n, m); mr = tuple(mr)
        if mr in oid: par[f(oid[m])] = f(oid[mr])
    lab = {m: keys[m] for m in lnas}
    bad = []
    for grp in qf.cospectralQuipuGroups(n).values():
        seeds = [tuple(int(c) for c in nk.QuipuAlgebra(*p).correspondingLNA().className()) for p in grp]
        name = [qf.formatQuipu(p) for p in grp]
        comp = {f(oid[s]): name[i] for i, s in enumerate(seeds)}
        assert len(comp) == len(seeds), 'seeds already joined'
        for m in lnas:
            if keys[m] == keys[seeds[0]]:
                c = f(oid[m])
                if c in comp: lab[m] = (keys[m], comp[c])
                else: bad.append(m); lab[m] = (keys[m], '?')
    return lab, bad
if __name__ == '__main__':
    n = int(sys.argv[1]); lab, bad = classes(n)
    print(n, 'LNAs', len(lab), 'classes', len(set(lab.values())), 'unresolved cospectral LNAs', len(bad), bad[:10])
