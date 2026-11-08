"""Round 121 maverick: is the n = 12 Coxeter key group of 47 LNAs (E-314; SNF(S) splits it 14 + 33) merged by the move census?
The census = `freeMoves.derivedOrbits` (rule table + free move + edge moves + double mutation), every move a recorded derived equivalence.
Restricted computation: every move preserves the Coxeter key, so an orbit of a group member never leaves the group; we therefore apply all
moves (rules, free strip, free ADD of a 2-relation = inverse strip, edges, doubles) to every one of the 47 members, union inside the group and
report (a) any move whose target lies OUTSIDE the key group (a key-changing "derived" move = a census bug), (b) the orbits, and (c) whether any
orbit contains both SNF classes (a merge across the split).  Also: a whole-length derivedOrbits cross-check of the same partition (CHECK=1).
SNF(S), S = C + C^T, computed twice: python-flint fmpz_mat.snf and pari matsnf.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/121/maverick_census47.py [CHECK=0|1]"""
import sys, os, collections, time
sys.path.insert(0, '.')
import flint, cypari
from quivermutation import coxeterTables as ct, freeMoves, lnaMoves, edgeMoves, doubleMutation
N = 12; TGT = (2, 2, 4, 0, 0, 3, 3, 0, 3, 0)
pari = cypari.pari


def snf_S(rl):
    C = ct.lnaCartanMatrix(N, tuple(rl)); S = [[C[i][j] + C[j][i] for j in range(N)] for i in range(N)]
    a = flint.fmpz_mat(S).snf(); d1 = tuple(sorted(abs(int(a[i, i])) for i in range(N)))
    d2 = tuple(sorted(abs(int(x)) for x in pari.matsnf(pari.matrix(N, N, [v for r in S for v in r]))))
    assert d1 == d2, (rl, d1, d2)
    return d1


def main():
    t0 = time.time(); keys = ct.lnaKeys(N)
    kk = keys[TGT]
    mem = sorted(tuple(r) for r, k in keys.items() if k == kk); memset = set(mem)
    print('group size', len(mem), 'key', kk, '%.0fs' % (time.time() - t0), flush=True)
    snf = {r: snf_S(r) for r in mem}; cls = collections.Counter(snf.values()); print('SNF(S) classes:', dict(cls), flush=True)
    lab = {r: ('A14' if cls[snf[r]] == 14 else 'B33') for r in mem}
    parent = {r: r for r in mem}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    edges = []; outside = []; kinds = collections.Counter()
    for r in mem:
        tg = [('rule', tuple(x)) for x in lnaMoves.rewritesOf(N, list(r))]
        tg += [('edge', tuple(x)) for x, _ in edgeMoves.rewritesOf(N, r)]
        tg += [('double', tuple(x)) for x, _ in doubleMutation.rewritesOf(N, r)]
        s = freeMoves.stripLengthTwo(r)
        if s != r: tg.append(('strip', s))
        for v in freeMoves.addableLengthTwo(N, r):
            d = list(r); d[v - 1] = 2; tg.append(('add2', tuple(d)))
        for kind, x in tg:
            kinds[kind] += 1
            if x not in memset:
                outside.append((r, kind, x, keys.get(x) == kk)); continue
            edges.append((r, kind, x)); parent[find(r)] = find(x)
    print('move applications by kind', dict(kinds), 'in-group edges', len(edges), 'targets outside group', len(outside), '%.0fs' % (time.time() - t0), flush=True)
    for o in outside[:10]: print('  OUTSIDE', o)
    cross = [(a, k, b) for a, k, b in edges if lab[a] != lab[b]]
    print('edges crossing the SNF split:', len(cross))
    for c in cross[:20]: print('  CROSS', c, snf[c[0]], snf[c[2]])
    orb = collections.defaultdict(list)
    for r in mem: orb[find(r)].append(r)
    print('orbits:', len(orb))
    for o in sorted(orb.values(), key=len, reverse=True):
        print('  size %2d labels %s' % (len(o), dict(collections.Counter(lab[r] for r in o))), ' '.join(ct.className(r) for r in o[:6]), '...' if len(o) > 6 else '')
    mixed = [o for o in orb.values() if len({lab[r] for r in o}) > 1]
    print('MIXED orbits (merge across the split):', len(mixed))
    # mirror (opposite algebra): not a derived equivalence, reported only to see whether the split is mirror-closed
    mir = collections.Counter((lab[r], lab.get(freeMoves.mirrorRow(N, r), 'out')) for r in mem); print('mirror label pairs', dict(mir))
    if os.environ.get('CHECK', '0') == '1':
        lnas, orbits = freeMoves.derivedOrbits(N, free=True, edges=True, doubles=True)
        full = [o for o in orbits.values() if set(o) & memset]
        print('whole-length derivedOrbits: orbits meeting the group', len(full), 'sizes', sorted(len(o) for o in full),
              'mixed', sum(1 for o in full if len({lab.get(r, 'out') for r in o}) > 1), '%.0fs' % (time.time() - t0), flush=True)
    print('done %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
