"""Round 029 (theorist): is the Coxeter key of each hand-built circuit algebra (D nn, W-type, G, H of rounds 023/027) a key of some LNA (or its dual) with the same number of vertices?
Usage: theorist_keys.py"""
import sys
sys.path.insert(0, '.'); MY = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec'))
sq = [(1,2),(1,3),(2,4),(3,4)]
cases = {
 'D nn (commute into both)': (sq+[(4,5),(4,6)], [[[1,2,4,5],[1,3,4,5]],[[1,2,4,6],[1,3,4,6]]]),
 'W-type (commute b1, kill b2 sides)': (sq+[(4,5),(4,6)], [[[1,2,4,5],[1,3,4,5]],[[1,2,4,6]],[[1,3,4,6]]]),
 'G': (sq+[(4,5),(1,7),(7,5)], [[[1,2,4,5],[1,7,5]],[[1,3,4,5],[1,7,5]]]),
 'H chain': ([(1,2),(1,3),(1,7),(2,4),(3,4),(7,4),(4,5),(4,6)], [[[1,2,4,5],[1,3,4,5]], [[1,3,4,6],[1,7,4,6]], [[1,2,4,6]], [[1,7,4,5]]]),
}
keysets = {}
def keys(n):
    if n not in keysets:
        s = set()
        for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
            for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
                s.add(search._coxeterKeyOrNone(alg))
        keysets[n] = s
    return keysets[n]
for name, (arr, rels) in cases.items():
    nodes = sorted({x for e in arr for x in e}); A = build(arr, rels, nodes)
    k = search._coxeterKeyOrNone(A)
    print('%-38s vertices %d key %s  in LNA keys at n=%d: %s' % (name, len(nodes), k, len(nodes), k in keys(len(nodes))))
