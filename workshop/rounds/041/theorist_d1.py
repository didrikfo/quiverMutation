"""Round 041 (theorist): check the closed form of the x^1 coefficient change for j = e_i, t = 0.
With G = C'^{-1}:  dP_1 := [x^1] det(x C_B + C_B^T) - [x^1] det(x C' + C'^T) = det(C') * ( G_vi - (G C'^T G)_{iv} - G_ii G_vv )   (indices: e_v row/col as in C_B = C' + E_{vi}).
Usage: theorist_d1.py dump.pkl   -> counts records with J = e_i, formula vs direct, and number with dP_1 = 0."""
import sys, pickle
import numpy as np
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']
xs = np.arange(0, n + 1); Vd = np.vander(xs.astype(float), n + 1, increasing=True)
def coef(C): return np.linalg.solve(Vd, np.array([np.linalg.det(t * C + C.T) for t in xs]))
ok = tot = zero = 0
for r in D['recs']:
    Jn = {i: j for i, j in r['J'].items() if j}
    if len(Jn) != 1 or list(Jn.values()) != [1]: continue
    V, v, CA, CB = r['V'], r['v'], r['CA'], r['CB']; k = V.index(v); i = V.index(list(Jn)[0])
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    Cp = rm @ CA @ rm.T; G = np.linalg.inv(Cp)
    direct = coef(CB)[1] - coef(Cp)[1]
    form = np.linalg.det(Cp) * (G[k, i] - (G @ Cp.T @ G)[i, k] - G[i, i] * G[k, k]) if abs(G[i, k]) < 1e-9 else None
    tot += 1
    if form is not None and abs(direct - form) < 1e-6: ok += 1
    if abs(direct) < 1e-6: zero += 1
print('records with J = e_i:', tot, 'formula agrees:', ok, 'dP_1 = 0:', zero)
