"""Round 042 (theorist): test the reduction on J != 0 records of a dump.
Also (b) F e_w = e_i and (c) F e_i = -e_m, F = Z Z^-T.
Facts tested per record (v = mutation vertex, w = its unique out-neighbour, i = the vertex with J_i = 1, Z = C_A on V\\v):
 (F1) out(v) = {w}; C'[a,v] = 0 for a != v (v last: C' = [[Z,0],[u^T,1]]); u = e_w - e_i;
 (F2) Q(x) = x[adjS_ii - adjS_wi - adjS_iw], S = xZ + Z^T, equals det(xC_B+C_B^T)-det(xC'+C'^T);
 (F3) low-order coefficients: Y_iw + Y_wi = 1 (Q_1 = 0), W_iw + W_wi - W_ii = 1 (Q_2 = 1), W = Y^T Z Y^T, Y = Z^-1.
Usage: theorist_reduce.py dump.pkl [onlyparentkey]"""
import sys, pickle
import numpy as np
import sympy as sp
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']; only = len(sys.argv) > 2
x = sp.symbols('x')
def poly(C): return sp.Poly(sp.Matrix(C.tolist()).applyfunc(lambda a: a).__mul__(1) * 0 + (x * sp.Matrix(C.tolist()) + sp.Matrix(C.T.tolist())), x) if False else sp.Poly((x * sp.Matrix(C.tolist()) + sp.Matrix(C.T.tolist())).det(), x)
cnt = Counter(); seenZ = set(); bad = []
for r in D['recs']:
    if only and not r.get("pkey", True): continue
    Jn = {i: j for i, j in r['J'].items() if j}
    if not Jn: continue
    V, v, CA, CB = r['V'], r['v'], r['CA'], r['CB']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    Cp = rm @ CA @ rm.T
    keep = [a for a in range(n) if a != k]
    cnt['J!=0 steps'] += 1; cnt['pkey', r.get('pkey')] += 1
    cnt['|supp J|', len(Jn)] += 1; cnt['dim J', tuple(sorted(Jn.values()))] += 1
    cnt['|out v|', len(r['out'])] += 1
    colv = all(Cp[a, k] == 0 for a in keep); cnt['C[a,v]=0', colv] += 1
    if len(r['out']) != 1 or len(Jn) != 1 or not colv: continue
    w = keep.index(V.index(r['out'][0])); i = keep.index(V.index(list(Jn)[0]))
    u = Cp[k, keep]; e = np.eye(n - 1, dtype=int)
    ok = bool((u == e[w] - e[i]).all()); cnt['u = e_w - e_i', ok] += 1
    if not ok: bad.append(r); continue
    Z = Cp[np.ix_(keep, keep)]
    key = (Z.tobytes(), w, i)
    if key in seenZ: cnt['duplicate (Z,w,i)'] += 1
    seenZ.add(key)
    S = x * sp.Matrix(Z.tolist()) + sp.Matrix(Z.T.tolist()); A = S.adjugate()
    Q = sp.expand(x * (A[i, i] - A[w, i] - A[i, w])); Q0 = sp.expand(poly(CB).as_expr() - poly(Cp).as_expr())
    cnt['F2 formula = Q', sp.expand(Q - Q0) == 0] += 1
    Pq = sp.Poly(Q, x); cs = Pq.all_coeffs()[::-1]; low = next(j for j, c in enumerate(cs) if c)
    cnt['lowest deg, coeff', (low, int(cs[low]))] += 1
    Y = np.rint(np.linalg.inv(Z)).astype(int); W = Y.T @ Z @ Y.T
    cnt['Y_iw+Y_wi', int(Y[i, w] + Y[w, i])] += 1; cnt['W_iw+W_wi-W_ii', int(W[i, w] + W[w, i] - W[i, i])] += 1
    F = Z @ Y.T; Fe_w = F[:, w]; Fe_i = F[:, i]
    cnt['(b) F e_w = e_i', bool((Fe_w == e[i]).all())] += 1
    cnt['(c) F e_i = -e_m', bool((-Fe_i >= 0).all() and (-Fe_i).sum() == 1)] += 1
    cnt['order: Y_wi=1 and Y_im=-1 (w<i<m)', bool(Y[w, i] == 1 and (-Fe_i).sum() == 1 and Y[i, int(np.argmax(-Fe_i))] == -1)] += 1
    Fi_ = np.rint(np.linalg.inv(F)).astype(int)
    def Pw(kk):
        M = np.eye(n - 1, dtype=int)
        for _ in range(abs(kk)): M = M @ (F if kk > 0 else Fi_)
        return M
    ss = tuple(sg for sg in range(-6, 7) if (Pw(sg)[:, w] == e[i]).all())
    cnt['s with F^s e_w = e_i (|s|<=6)', ss] += 1
    cc = {kk: int((Y.T @ Pw(kk))[w, w]) for kk in range(-8, 9)}
    cnt['Serre c_-k = c_k-1 (k<=6)', all(cc[-kk] == cc[kk - 1] for kk in range(0, 7))] += 1
    if ss:
        sg = ss[0]; rho = lambda kk: cc[kk] - cc[kk + sg] - cc[kk - sg]
        cnt['moment formula: rho_0 = 0 and -rho_1 = 1', (rho(0), -rho(1))] += 1
    cnt['Z_iw, Z_wi', (int(Z[i, w]), int(Z[w, i]))] += 1; cnt['Y_iw, Y_wi', (int(Y[i, w]), int(Y[w, i]))] += 1
for kx in sorted(cnt, key=str): print(kx, cnt[kx])
print('records with u != e_w - e_i:', len(bad))
