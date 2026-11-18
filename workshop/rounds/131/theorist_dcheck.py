"""Round 131, theorist, answer to referee item 2 and 3.

(d) per row: replays `theorist_thin.py sweep N seed` (same RNG stream, mixed
sample only) and at every gate-admitted vertex with J != 0 checks
  (d1) some i != v with d_i >= 2 and some out-arrow beta such that
       e_i I e_t(beta) is not spanned by paths;
  (d2) the same, with i restricted to the i where J_i != 0.
e_i I e_j is spanned by paths iff its fully reduced echelon basis consists of
single paths (a coordinate subspace has a unique RREF of unit vectors).
Control: the same test at admitted J = 0 vertices with d >= 2 (informative).

bag: E-329's counterexample to the library gate (toolsmith_cycgate.py bag),
recomputed with the independent engine: d_w per w at v = 1 and the gate.

usage: theorist_dcheck.py rows [N] [seed] | bag
"""
import random, sys, warnings
from theorist_thin import Alg, randomAlg, paths

def pathSpanned(alg, i, j):
    return all(len(r) == 1 for r in alg.pivots(i, j).values())

def rows(N, seed):
    rng = random.Random(seed)
    tot = d1 = d2 = 0; ctl = ctl_nm = 0; bad = []
    for k in range(N):
        alg = randomAlg(rng, False)
        for v in range(alg.n):
            if not alg.out(v): continue
            g = alg.gate(v)
            ds = {i: alg.dJ(i, v) for i in range(alg.n) if i != v}
            dmax = max(d for d, _ in ds.values()); J = sum(j for _, j in ds.values())
            if not g or dmax < 2: continue
            nm = lambda i: ds[i][0] >= 2 and any(not pathSpanned(alg, i, b[1]) for b in alg.out(v))
            if J > 0:
                tot += 1
                a = any(nm(i) for i in ds); b = any(nm(i) for i in ds if ds[i][1] > 0)
                d1 += a; d2 += b
                if not b: bad.append((k, v))
            else:
                ctl += 1; ctl_nm += any(nm(i) for i in ds)
    print(f"(2,T,T) rows {tot}: (d1) some i with d_i>=2 and non-monomial e_i I e_t(beta): {d1}; "
          f"(d2) at an i with J_i != 0: {d2}; failures {bad}")
    print(f"control (2,T,F) rows {ctl}: non-monomial piece present at {ctl_nm} (informative, not claimed either way)")

def bag():
    # E-329 / toolsmith_cycgate.py bag, vertices 1..5 -> 0..4:
    # a:4>5, b:5>4, g:5>1, x1:1>2, x2:2>3; relations ab, bagx1; v = 1
    a, b, g, x1, x2 = (3, 4, 'a'), (4, 3, 'b'), (4, 0, 'g'), (0, 1, 'x1'), (1, 2, 'x2')
    alg = Alg(5, [a, b, g, x1, x2], [{(a, b): 1}, {(b, a, g, x1): 1}])
    v = 0
    for w in range(5):
        if w == v: continue
        std = alg.std(w, v)
        print(f"w={w+1}: d_w={len(std)} basis={[''.join(x[2] for x in p) for p in std]}")
    # gate with the full (non-simple) quantifier
    ref = [p for i in range(5) if i != v for p in alg.std(i, v)
           if not any(alg.nf({p + (bb,): 1}) for bb in alg.out(v))]
    simple = lambda p: len({x[0] for x in p} | {p[-1][1]}) == len(p) + 1
    print("refusing nonzero paths:", [(''.join(x[2] for x in p), 'simple' if simple(p) else 'non-simple') for p in ref],
          "; exact gate:", alg.gate(v))

if __name__ == "__main__":
    warnings.simplefilter("ignore")
    if sys.argv[1] == "rows":
        rows(int(sys.argv[2]) if len(sys.argv) > 2 else 600, int(sys.argv[3]) if len(sys.argv) > 3 else 131)
    else:
        bag()
