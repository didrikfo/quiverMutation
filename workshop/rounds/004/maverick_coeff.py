"""Is the second Coxeter-polynomial coefficient (T^(n-2)) a function of (cords, relations) or of anything simple?
Over every quipu quiver with monomial relations of order N (all orientations, all ideals, no lines)."""
import sys, collections, numpy as np
from quivermutation import quipuRelations as qr
n = int(sys.argv[1]); minA = 2
tab = collections.defaultdict(collections.Counter)
for par, edges, ori, aut in qr.orientedQuipus(n):
    if qr.isLinearlyOriented(n, edges, ori): continue
    succ = qr.successors(n, edges, ori); paths = qr.directedPaths(n, succ)
    base, cand, kill, comp = qr.relationData(n, paths, minA)
    cords = sum(1 for x in par[1] if x > 0)
    for mask, chosen in qr.relationSets(base, kill, comp):
        C = np.array(qr.cartanFromMask(mask, n), dtype=float)
        # Phi = -C^-T C ; coefficients of char poly
        Phi = -np.linalg.inv(C).T @ C
        cp = np.rint(np.poly(Phi)).astype(int)     # leading 1, then -tr, then e2 ...
        tab[(cords, len(chosen))][(cp[1], cp[2])] += 1
for k in sorted(tab): print("cords,rels", k, "(c_{n-1}, c_{n-2}) values:", dict(tab[k]))
