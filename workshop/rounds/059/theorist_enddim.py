"""Round 059 (theorist), T10/T5: dim-level test of the lemma  End_K(T_v) = H_0 of the mutated dg quiver = repo rewrite, at EVERY vertex, J = 0 or not.
For every LNA and dual of length N and every vertex v with an arrow out (no gate): child = procedure.mutateAtVertex (steps 1-7, unreduced);
H[i][j] = dim Hom_K(T_i,T_j) from the independent complex code of round 050; compare labelled Cartan(child) with H (and H^T), and with R C R^T.
Prediction (lemma): Cartan(child) == H or H^T at all steps; Cartan(child) - R C R^T = sum_{m>=1} (-1)^{m+1} dim Hom(T,T[-m]) entrywise, so = 0 iff
(at the level of dimensions) no negative-degree Homs; J != 0 steps differ.
Usage: .venv/bin/python workshop/rounds/059/theorist_enddim.py N"""
import sys
sys.path.insert(0, '.')
from collections import Counter
_st = open('workshop/rounds/050/skeptic_tilt.py').read()
exec(_st[_st.index('from fractions'):_st.index('def stepTest')])
exec(_st[_st.index('def stepTest'):_st.index('def stepTest') + 400].split('\n\nclass')[0])
from quivermutation import nakayama as nk, pathAlgebra
N = int(sys.argv[1])
cnt = Counter(); bad = []
for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
    for tag, alg in (('lna', lna), ('dual', pathAlgebra.dualPathAlgebra(lna))):
        q = alg.quiver; rels = procedure.relationsFrom(alg)
        V = sorted(q.nodes)
        for v in V:
            if not ap.arrowsOutOf(q, v): continue
            try:
                nq, nr = procedure.mutateAtVertex(q, rels, v)
                C = ap.cartanMatrix(nq, nr, exact=True)
                D = procedure.cartanDiscrepancy(q, rels, v, nq, nr)
            except Exception as e:
                cnt['error'] += 1; continue
            _, res = stepTest(alg, v, -1)
            H, Hm, Hp = res[0], res[-1], res[1]
            n = len(V)
            eqH = all(C[i][j] == H[i][j] for i in range(n) for j in range(n))
            eqHt = all(C[i][j] == H[j][i] for i in range(n) for j in range(n))
            J0 = sum(map(sum, Hm)) == 0 and sum(map(sum, Hp)) == 0
            # Euler: R C R^T - Cartan(child) should equal H^{-1}-type correction; D = RCR^T - C
            eul = all(D[i][j] == 0 for i in range(n) for j in range(n))
            # D[i][j] = (RCR^T)_{ij} - C_ij = chi - h0 = -(h^{-1} - h^{-2} ...) ; test vs -Hm (valid when h^{-2}=0)
            dmH = all(D[i][j] == -Hm[i][j] for i in range(n) for j in range(n)); dmHt = all(D[i][j] == -Hm[j][i] for i in range(n) for j in range(n))
            dm = dmH or dmHt
            cnt[('SPLIT', 'D=-Hm: ' + ('H' if dmH else '') + ('Ht' if dmHt else '') , 'Hp==0' if sum(map(sum, Hp)) == 0 else 'Hp!=0', 'Hm==0' if sum(map(sum, Hm)) == 0 else 'Hm!=0', 'eqH' if eqH else '', 'eqHt' if eqHt else '')] += 1
            cnt[(tag, 'J0' if J0 else 'Jne0', ('End=child(H)' if eqH else 'End=child(H^T only)') if (eqH or eqHt) else 'MISMATCH', 'euler0' if eul else 'euler!=0', 'D=-Hm' if dm else 'D!=-Hm')] += 1
            if not (eqH or eqHt): bad.append((lna, tag, v))
for k, c in sorted(cnt.items(), key=str): print(c, k)
print('mismatches', len(bad), bad[:5])
