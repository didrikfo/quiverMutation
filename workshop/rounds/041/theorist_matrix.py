"""Round 041 (theorist): pure matrix-level search, no algebra. Does C'+E_{vi} (or E_{iv}) ever have det(xC+C^T)/det C equal to that of C'
for C' = r C r^T, C unitriangular (entries 0..m), r = identity except row v = -e_v + sum_{w in S} e_w (S a subset of the other vertices, multiplicity 1)?
Also drops the requirement that C' come from r: arbitrary integer unimodular C' with the same row-v shape is NOT tried.
Mode CON=2 additionally takes S = the 'minimal' support {w: C[v,w]>=1, no u with C[v,u]>=1 and C[u,w]>=1} (arrows if no relation hides a path) and demands sum_{w in S} C[i,w] >= C[i,v]-1 (rank bound for J_i = 1).
Mode CON=1 (env): realistic constraints: S subset of {w: C[v,w]>=1}, only the perturbation E_{v i} with C[i,v]>=2 (d_i >= 2, J_i = 1 <= d_i - 1).
Usage: [CON=1] theorist_matrix.py n m [samples]   (samples = 0: exhaustive)"""
import sys, itertools, random, os
CON = os.environ.get('CON') in ('1', '2'); CON2 = os.environ.get('CON') == '2'
import numpy as np
n, m = int(sys.argv[1]), int(sys.argv[2]); samples = int(sys.argv[3]) if len(sys.argv) > 3 else 0
xs = np.arange(0, n + 1); Vd = np.vander(xs.astype(float), n + 1, increasing=True)
def poly(C): return tuple(int(round(c)) for c in np.linalg.solve(Vd, np.array([np.linalg.det(t * C + C.T) for t in xs])))
pos = [(a, b) for a in range(n) for b in range(a + 1, n)]
def mats():
    if samples:
        for _ in range(samples): yield {p: random.randint(0, m) for p in pos}
    else:
        for vals in itertools.product(range(m + 1), repeat=len(pos)): yield dict(zip(pos, vals))
tot = hit = 0; ex = []
for d in mats():
    C = np.eye(n, dtype=int)
    for (a, b), c in d.items(): C[a, b] = c
    for v in range(n):
        if CON2:
            sup = [w for w in range(n) if w != v and C[v, w] >= 1]
            Smin = tuple(w for w in sup if not any(u != w and C[u, w] >= 1 for u in sup))
            Slist = [Smin]
        else:
            Slist = itertools.chain.from_iterable(itertools.combinations([w for w in range(n) if w != v and (not CON or C[v, w] >= 1)], s) for s in range(0, n))
        for S in Slist:
            r = np.eye(n, dtype=int); r[v, v] = -1
            for w in S: r[v, w] = 1
            Cp = r @ C @ r.T; P0 = poly(Cp)
            for i in range(n):
                if i == v: continue
                for (a, b) in (((v, i),) if CON else ((v, i), (i, v))):
                    if CON and C[i, v] < 2: continue
                    if CON2 and sum(C[i, w] for w in S) < C[i, v] - 1: continue
                    E = np.zeros((n, n), dtype=int); E[a, b] = 1; tot += 1
                    if poly(Cp + E) == P0:
                        hit += 1
                        if len(ex) < 5: ex.append((C.tolist(), v, S, (a, b)))
print('n', n, 'm', m, 'trials', tot, 'preserving', hit)
for e in ex: print(e)
