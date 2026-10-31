"""Round 103 literature check: arXiv:2609.39536 (Deng-Geng-Liu) Theorem 'conditions for approximation is injective'
(item (3): for every outgoing beta at k there is an incoming alpha at k with alpha*beta not in I; here alpha*beta = alpha then beta)
against the repo's gate `procedure.isMutable`, on every gentle LNA (all relations of length 2: row entries in {0,2}), n = 3..N.
Row convention: n-2 entries, entry j-1 = relation of j -> j+1 -> j+2 starting at vertex j (length 2).
Paper criterion A(v) = v == n or (v >= 2 and r_{v-1} == 0);  mirror B(v) = v <= n-1 and (v == 1 or r_{v-1} == 0).
Compare A with the repo LEFT gate and B with the repo RIGHT gate (and the crossed pairs, to see the control).
usage: timeout 10m .venv/bin/python workshop/rounds/103/experimentalist_gentle_gate.py [N=9]"""
import sys, itertools
from quivermutation import nakayama, procedure, lnaMoves as lm
N = int(sys.argv[1]) if len(sys.argv) > 1 else 9
tot = {}
for n in range(3, N + 1):
    cnt = {'rowsxverts': 0, 'R~B': 0, 'L~A': 0, 'R~A': 0, 'L~B': 0, 'R_pass': 0, 'L_pass': 0}
    for r in itertools.product([0, 2], repeat=max(n - 2, 0)):
        alg = nakayama.LinearNakayamaAlgebra(n, list(r))
        q, rels = alg.quiver, procedure.relationsFrom(alg)
        oq, orels = procedure._opposite(q, rels)
        for v in range(1, n + 1):
            R = bool(lm._quiet(procedure.isMutable, q, rels, v))
            L = bool(lm._quiet(procedure.isMutable, oq, orels, v))
            rv1 = r[v - 2] if 2 <= v <= n - 1 else 0     # r_{v-1}; vertex v-1 must be in 1..n-2
            A = (v == n) or (v >= 2 and rv1 == 0)
            B = (v <= n - 1) and (v == 1 or rv1 == 0)
            cnt['rowsxverts'] += 1
            cnt['R~B'] += R == B; cnt['L~A'] += L == A; cnt['R~A'] += R == A; cnt['L~B'] += L == B
            cnt['R_pass'] += R; cnt['L_pass'] += L
    print(n, cnt, flush=True)
