"""Round 054 (theorist), T10: is generation of K^b(proj A) by T = tiltingPlus complex automatic?
Argument (see theorist.md): T_v = cone(P_v -> +_{v->h} P_h); acyclic => every h != v, so M=+P_h lies in add(T) and P_v = cone(M -> T_v)[-1] in thick(T);
all P_i, i != v are summands; hence thick(T) = K^b(proj A).  The script checks the HYPOTHESES of that argument and the m=+1 / m=-1 split on every step:
  (H1) quiver acyclic (no loop, no 2-cycle), (H2) arrows out of v = all of ap.arrowsOutOf, targets h != v and distinct summands P_h in T,
  (H3) T has n pairwise non-isomorphic indecomposable summands (T_v not iso to P_i: its terms differ in degree), (H4) class matrix X of T in K0 has det +-1,
  (H5) the Euler form of T: sum_m (-1)^m dim Hom(T_i,T_j[m]) = X C X^T (consistency, computed from the same Hom code),
  and records dim Hom(T,T[1]) and dim Hom(T,T[-1]) (J = 0 is the vanishing of both).
Usage:  .venv/bin/python workshop/rounds/054/theorist_gen.py lna N      all LNAs and duals of length N, every vertex with mutationIsPossibleAtVertex
        .venv/bin/python workshop/rounds/054/theorist_gen.py path13     the 13 E-163 edges (needs /tmp/tsm/c1.pkl, see toolsmith_endt_run.py)"""
import sys
sys.path.insert(0, '.')
MODE = sys.argv[1] if len(sys.argv) > 1 else 'lna'
ARGS = sys.argv[2:]
from fractions import Fraction
import networkx as nx
from quivermutation import arrowPaths as ap, procedure, mutation, reduction, pathAlgebra
from quivermutation import nakayama as nk
_st = open('workshop/rounds/050/skeptic_tilt.py').read()
P = 2147483647
exec(_st[_st.index('class Alg'):_st.index('def stepTest')])

def check(alg, v):
    A = Alg(alg); V = sorted(alg.quiver.nodes); n = len(V)
    out = {}
    out['H1_acyclic'] = not any(True for _ in nx.simple_cycles(alg.quiver))
    outs = ap.arrowsOutOf(A.q, v)
    out['H2_targets_ne_v'] = all(b[1] != v for b in outs) and len(outs) > 0
    T = complexT(A, v, V, -1)
    # summands of T: P_i (i != v) in degree 0 and T_v with terms P_v (deg -1), P_h (deg 0). Non-iso: T_v has a degree -1 term.
    out['H3_n_summands'] = (len(T) == n) and all(set(T[i][0]) == {0} for i in V if i != v) and (-1 in T[v][0])
    out['H2b_P_h_in_addT'] = all(h in V and h != v for h in [b[1] for b in outs])
    # K0 class matrix: [T_i] = e_i, [T_v] = sum_h e_h - e_v  -> det
    import sympy
    X = sympy.zeros(n, n)
    for r, i in enumerate(V):
        if i != v: X[r, V.index(i)] = 1
        else:
            for b in outs: X[r, V.index(b[1])] += 1
            X[r, V.index(v)] -= 1
    out['H4_det'] = int(X.det())
    # Euler form from Hom dims
    hom = {m: [[homdim(A, T[i], T[j], m) for j in V] for i in V] for m in (-1, 0, 1)}
    out['hom_p1'] = sum(map(sum, hom[1])); out['hom_m1'] = sum(map(sum, hom[-1]))
    return out

def summarize(rows):
    keys = ['H1_acyclic', 'H2_targets_ne_v', 'H2b_P_h_in_addT', 'H3_n_summands']
    print('steps', len(rows), {k: sum(1 for r in rows if r[k]) for k in keys}, 'det=+-1:', sum(1 for r in rows if abs(r['H4_det']) == 1))
    from collections import Counter
    print('(hom_p1==0, hom_m1==0) counts:', Counter((r['hom_p1'] == 0, r['hom_m1'] == 0) for r in rows))

if MODE == 'lna':
    N = int(ARGS[0]); rows = []
    for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
            for v in sorted(alg.quiver.nodes):
                if not ap.arrowsOutOf(alg.quiver, v): continue
                rows.append(check(alg, v))
    summarize(rows)
elif MODE == 'path13':
    sys.argv = ['x', '/tmp/tsm/c1.pkl', '1', '5', '6', '400000', 'paths']
    import pickle
    src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
    exec(compile(src, 'tp', 'exec'))
    recs = [r for r in pickle.load(open('/tmp/tsm/c1.pkl', 'rb')) if r['kind'] == 'fail']
    rows = []; x0 = recs[14]['childObj']
    for name, x, mv in [('child', x0, 'F1 F3 F1 F5 F4 R7 R1'.split()), ('lna9', classes[base][9], 'F7 R2 R1 F2 F7 R3'.split())]:
        for s_ in mv:
            kind, v = s_[0], int(s_[1:]); a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
            r = check(a, v); rows.append(r); print(name, s_, r, flush=True)
            c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v)); x = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
    summarize(rows)
    # failing steps (J != 0): which sign fails?
    frows = [check(r['parentObj'], r['v']) for r in recs]
    print('FAILING steps (16):'); summarize(frows)
    print([(f['hom_p1'], f['hom_m1']) for f in frows])
