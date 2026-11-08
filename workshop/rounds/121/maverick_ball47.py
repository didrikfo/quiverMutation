"""Round 121 maverick: a key-guarded mutation ball from one LNA of the 14-class of the n = 12 key group of 47 (E-314), as a T10 / H-015 instrument.
Walk: the library's own `search.mutationSearchDepthFirst` (gate + Coxeter-key guard, deduped by `fingerprint.Visited`), from the root and from its
opposite (as `linesReachedFrom(alsoDual=True)` does).  Detector: at EVERY node reached (not only LNAs), SNF(C + C^T) of the integer Cartan matrix.
SNF of the symmetrised Euler form is a derived invariant for acyclic algebras (finite gl.dim; Ladkani chain, E-308), and it is (0,1^10,20) on the
14-class and (0,1^10,12) on the 33-class.  A node with SNF != root's is a node the guarded walk reached OUTSIDE the root's derived class (H-015 NO).
Also recorded: the LNAs reached (line nodes), labelled A14 / B33 / other.
Controls (CTRL env):  'planted' -- the detector is shown the 33-class LNAs as if they were ball nodes and must flag every one (planted crossing);
                      'off'     -- the same ball with coxeterGuard=False: nodes whose key moved are certainly NOT derived equivalent to the root;
                                   reports how many of them the SNF detector flags, and how many key-EQUAL nodes with a different SNF the unguarded walk meets.
usage (repo root): ROOT=3040500030 DEPTH=3 [SIDES=right,opposite] [CTRL=off|planted] timeout 10m .venv/bin/python -W ignore workshop/rounds/121/maverick_ball47.py"""
import sys, os, copy, time, collections
sys.path.insert(0, '.')
import flint
import networkx as nx
from quivermutation import coxeterTables as ct, search, fingerprint, invariants, pathAlgebra, nakayama as nk
N = 12; TGT = (2, 2, 4, 0, 0, 3, 3, 0, 3, 0)
ROOT = os.environ.get('ROOT', '3040500030'); DEPTH = int(os.environ.get('DEPTH', '3')); CTRL = os.environ.get('CTRL', '')
_snf = {}


def snf_of(C):
    k = tuple(map(tuple, C))
    if k not in _snf:
        n = len(C); S = flint.fmpz_mat([[C[i][j] + C[j][i] for j in range(n)] for i in range(n)]).snf()
        _snf[k] = tuple(sorted(abs(int(S[i, i])) for i in range(n)))
    return _snf[k]


def line_rl(alg):
    q = alg.quiver
    if q.number_of_edges() != N - 1 or len(q.nodes) != N or not nx.is_directed_acyclic_graph(q) or nx.dag_longest_path_length(q) != N - 1: return None
    L = search._renumberedLine(alg); rl = [0] * (N - 2)
    for rel in L.rels:
        if len(rel) != 1: return None
        p = rel[0]; rl[p[0] - 1] = len(p) - 1
    return tuple(rl)


def main():
    t0 = time.time(); keys = ct.lnaKeys(N); kk = keys[TGT]
    mem = [tuple(r) for r, k in keys.items() if k == kk]
    snfL = {r: snf_of(ct.lnaCartanMatrix(N, r)) for r in mem}; cnt = collections.Counter(snfL.values())
    lab = {r: ('A14' if cnt[snfL[r]] == 14 else 'B33') for r in mem}
    root_rl = tuple(int(c) for c in ROOT); assert root_rl in lab, 'root not in the 47-group'
    base = snfL[root_rl]; print('root', ROOT, lab[root_rl], 'SNF tail', base[-1], 'depth', DEPTH, 'ctrl', CTRL or 'none', flush=True)
    if CTRL == 'planted':
        flagged = sum(snf_of(invariants.integerCartanMatrix(nk.LinearNakayamaAlgebra.fromClassName(ct.className(r)))) != base for r in mem if lab[r] == 'B33')
        same = sum(snf_of(invariants.integerCartanMatrix(nk.LinearNakayamaAlgebra.fromClassName(ct.className(r)))) != base for r in mem if lab[r] == 'A14')
        print('PLANTED: 33-class LNAs flagged %d of 33; 14-class LNAs flagged %d of 14 (want 33, 0)' % (flagged, same)); return
    guard = CTRL != 'off'; stats = collections.Counter(); flaggedN = []; lnas = collections.Counter(); snfhist = collections.Counter()
    for side in os.environ.get('SIDES', 'right,opposite').split(','):
        start = nk.LinearNakayamaAlgebra.fromClassName(ROOT)
        start = copy.deepcopy(start) if side == 'right' else pathAlgebra.dualPathAlgebra(start)
        vis = fingerprint.Visited()

        def visit(alg, path, side=side):
            stats['visits'] += 1
            k = search._coxeterKeyOrNone(alg)
            if k is None: stats['cyclic/no key'] += 1; return
            C = invariants.integerCartanMatrix(alg); s = snf_of(C); snfhist[(k == kk, s[-1] if s[0] == 0 else ('nonsing', s[-1]))] += 1
            if s != base:
                stats['SNF differs' + (' (key moved)' if k != kk else ' (KEY KEPT)')] += 1
                if k == kk and len(flaggedN) < 20: flaggedN.append((side, list(path), s))
            elif k != kk: stats['key moved, SNF same'] += 1
            rl = line_rl(alg)
            if rl is not None:
                if side == 'opposite': rl = tuple(rl)  # opposite LNA, labelled as found (the 47-group is mirror-closed, label preserved)
                lnas[(ct.className(rl), lab.get(rl, 'other'))] += 1
        search.mutationSearchDepthFirst(start, DEPTH, printOutput=False, visitor=visit, coxeterGuard=guard, visited=vis)
        print(' side', side, 'visits so far', stats['visits'], 'distinct expanded', len(vis.depths), 'snf cache', len(_snf), '%.0fs' % (time.time() - t0), flush=True)
    print('stats', dict(stats))
    print('SNF histogram (key kept?, SNF tail):', dict(snfhist))
    print('LNAs reached:', len(lnas), 'labels', dict(collections.Counter(l for _, l in lnas)))
    for (name, l), c in sorted(lnas.items()): print('  ', name, l, c)
    for f in flaggedN: print('  KEY-KEPT SNF CHANGE', f)
    print('VERDICT: key-kept nodes outside the root derived class:', stats['SNF differs (KEY KEPT)'], '%.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
