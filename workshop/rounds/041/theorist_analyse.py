"""Round 041 (theorist): matrix-level analysis of dumped steps. Usage: theorist_analyse.py dump.pkl
For each record: C' = r C_A r^T, check C_B == C' + H (H_{v i} = J_i); polynomial P_B(x)=det(xC_B+C_B^T) vs P'(x); classify J support;
for J = e_i (single support, dim 1): D(x) = P_B - P' ; also hypothetical test: for every real (C_A, v) and every i != v, does C'+E_{vi} (or E_{iv}) have the same polynomial as C'?"""
import sys, pickle
import numpy as np
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); recs = D['recs']; n = D['n']
def poly(C):
    # det(xC+C^T) is a degree-n integer polynomial; return its exact coefficients by interpolation at x = 0..n (floats rounded)
    xs = np.arange(0, n + 1); ys = [round(np.linalg.det(t * C + C.T)) for t in xs]
    V = np.vander(xs.astype(float), n + 1, increasing=True)
    return tuple(int(round(c)) for c in np.linalg.solve(V, np.array(ys, dtype=float)))
cong_bad = 0; sup = Counter(); seenCp = {}; hyp = Counter(); rows = []
for r in recs:
    CA, CB, V, v = r['CA'], r['CB'], r['V'], r['v']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    Cp = rm @ CA @ rm.T
    H = np.zeros((n, n), dtype=int)
    for i, jj in r['J'].items(): H[k, V.index(i)] = jj
    if not (CB == Cp + H).all(): cong_bad += 1
    Jn = {i: jj for i, jj in r['J'].items() if jj}
    sup[(len(Jn), tuple(sorted(Jn.values())), r['passed'])] += 1
    key = (Cp.tobytes(), k)
    if key not in seenCp:
        seenCp[key] = (Cp, k)
print('records', len(recs), 'C_B != C\'+H:', cong_bad)
print('(support size, dims, passed): count'); [print('  ', a, b) for a, b in sorted(sup.items(), key=str)]
print('distinct (C\', v):', len(seenCp))
# hypothetical single-entry perturbations
base = {}
for (kb, (Cp, k)) in seenCp.items():
    P0 = poly(Cp)
    for i in range(n):
        if i == k: continue
        E = np.zeros((n, n), dtype=int); E[k, i] = 1
        P1 = poly(Cp + E)
        hyp[('row-v entry (v,i)', P1 == P0)] += 1
        E2 = np.zeros((n, n), dtype=int); E2[i, k] = 1
        hyp[('entry (i,v)', poly(Cp + E2) == P0)] += 1
print('hypothetical C\'+E_{vi} / C\'+E_{iv} with same Coxeter-type polynomial as C\':'); [print('  ', a, b) for a, b in sorted(hyp.items(), key=str)]
