"""Round 031 (theorist): hand-built algebras whose J_i (kernel of g_i at v) has a circuit of length 3 (T3) or an nn 2-cycle (D), plus pendant
vertices added anywhere; is the Coxeter key in the LNA/dual-LNA key set of the same length?  Also verifies the identity J_i = socle part
{x in e_iAe_v : x*rad = 0}  (the kernel is the module ker(P_v -> sum P_tb) = Hom(S_v, e_iA)).
Usage: theorist_t3.py [maxextra]"""
import sys, itertools
sys.path.insert(0, '.'); MY = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec'))
maxextra = int(MY[1]) if len(MY) > 1 else 2
keysets = {}
def keys(n):
    if n not in keysets:
        s = set()
        for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
            for alg in (lna, pathAlgebra.dualPathAlgebra(lna)): s.add(search._coxeterKeyOrNone(alg))
        keysets[n] = s
    return keysets[n]
def probe(arr, rels, v):
    nodes = sorted({x for e in arr for x in e}); A = build(arr, rels, nodes)
    kd, mono = kerdim(A, v, procedure.relationsFrom(A))
    return A, kd, mono, search._coxeterKeyOrNone(A), len(nodes)
# i=1, a1,a2,a3 = 2,3,4, v=5, t1=6, t2=7
tri = [(1,2),(1,3),(1,4),(2,5),(3,5),(4,5),(5,6),(5,7)]
T3 = [[[1,2,5,6]], [[1,4,5,7]], [[1,3,5,7],[1,2,5,7]], [[1,3,5,6],[1,4,5,6]]]
sq = [(1,2),(1,3),(2,5),(3,5),(5,6),(5,7)]
D = [[[1,2,5,6],[1,3,5,6]], [[1,2,5,7],[1,3,5,7]]]
W = [[[1,2,5,6],[1,3,5,6]], [[1,2,5,7]], [[1,3,5,7]]]
base = {'T3 (circuit 3, ground-ground path)': (tri, T3), 'D (nn)': (sq, D), 'W (control)': (sq, W)}
tot = {}
for name, (arr, rels) in base.items():
    A, kd, mono, k, n = probe(arr, rels, 5)
    print('%-40s n=%d kerdim=%d single-path-in-J=%s key in LNA keys: %s' % (name, n, kd, mono, k in keys(n)), flush=True)
    # add up to maxextra pendant vertices (arrow to/from an existing or new vertex), no new relations, never an arrow out of v or into i-side creating new i->v paths
    nodes0 = sorted({x for e in arr for x in e}); hits = 0; tried = 0
    def ext(arr_, nodes_, depth):
        global hits, tried
        if depth == maxextra: return
        new = max(nodes_) + 1
        for u in nodes_:
            for e in ((u, new), (new, u)):
                if e[1] == 1 and False: continue
                if e[0] == 5 and e[1] != 6 and e[1] != 7: pass
                arr2 = arr_ + [e]; nodes2 = nodes_ + [new]
                if e[0] == 5 or (e[1] == 5) or (e[1] == 1 and False): continue  # keep out-arrows of v and paths into v unchanged
                B, kd2, mono2, k2, n2 = probe(arr2, rels, 5)
                tried += 1
                if kd2 > 0 and not mono2:
                    if k2 in keys(n2): hits += 1; print('   LNA KEY HIT', arr2, 'n =', n2, flush=True)
                ext(arr2, nodes2, depth + 1)
    ext(arr, nodes0, 0)
    print('   extensions tried %d, LNA-key hits with J != 0 and no single path: %d' % (tried, hits), flush=True)
