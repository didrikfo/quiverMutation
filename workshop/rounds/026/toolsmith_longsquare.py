"""Round 026 (toolsmith): `longSquare` that handles parallel arrows (E-108).

`longSquareOld(alg, v)` is the round 023/022 test (scholar_longsquare.py), on the vertex-list `alg.rels`: two paths through a
doubled arrow project to the same vertex list, so `len({q[-3]}) == len(rel)` fails and a real long square is missed.
`longSquare(alg, v)` reads the arrow relations (`procedure.relationsFrom`, which prefers `arrowRels`) and asks for DISTINCT
second-to-last ARROWS instead of distinct predecessor vertices. Same meaning otherwise: out-degree of v is 1, arrow b = v -> e;
some relation has >= 2 terms, every path has >= 2 arrows, ends with b, and all start at one vertex. On quivers without parallel
arrows the two agree (checked: `.venv/bin/python workshop/rounds/026/toolsmith_longsquare.py`, n = 7 one step out of every LNA).
Import: importlib from this path (tests/test_longsquare.py does)."""
import sys
sys.path.insert(0, '.')
from quivermutation import arrowPaths as ap, procedure


def longSquareOld(alg, v):
    outs = ap.arrowsOutOf(alg.quiver, v)
    if len(outs) != 1: return False
    e = outs[0][1]
    for rel in alg.rels:
        if len(rel) >= 2 and all(len(q) >= 3 and q[-1] == e and q[-2] == v for q in rel) and len({q[0] for q in rel}) == 1 and len({q[-3] for q in rel}) == len(rel):
            return True
    return False


def longSquare(alg, v):
    outs = ap.arrowsOutOf(alg.quiver, v)
    if len(outs) != 1: return False
    b = outs[0]
    for rel in procedure.relationsFrom(alg):
        paths = list(rel)
        if len(paths) >= 2 and all(len(q) >= 2 and q[-1] == b for q in paths) \
                and len({q[0][0] for q in paths}) == 1 and len({q[-2] for q in paths}) == len(paths):
            return True
    return False


if __name__ == '__main__':
    from quivermutation import nakayama as nk, mutation, reduction, pathAlgebra
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    tot = diff = pos = 0; par = 0
    for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
            cands = [(alg, v) for v in alg.vertices()]
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, v))
                cands += [(ch, w) for w in ch.vertices()]
            for A, w in cands:
                tot += 1; o, nw = longSquareOld(A, w), longSquare(A, w); pos += nw
                if o != nw:
                    diff += 1; par += A.hasParallelArrows()
    print('n', n, 'tests', tot, 'new True', pos, 'old != new', diff, '(of which on a parallel-arrow algebra:', par, ')')
