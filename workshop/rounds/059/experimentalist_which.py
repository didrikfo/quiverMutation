"""E-153 follow-up: for every LNA of length n record the unique right v (and left v) whose
mutation gives an LNA, and the new relation vector.  Usage: which.py N [--print]  -> table + JSON-ish lines."""
import sys, collections, json
from quivermutation import lnaMoves as lm, mutation, nakayama, pathAlgebra
n = int(sys.argv[1]); show = '--print' in sys.argv

def allLNAs(n):
    out = []
    def rec(j, cur):
        if j == n - 2:
            if lm.isAdmissible(n, cur): out.append(list(cur))
            return
        for a in [0] + list(range(2, n - 1 - j + 1)):
            cur.append(a); rels = lm.relationsOf(cur)
            if all(e[0] + e[1] < l[0] + l[1] for e, l in zip(rels, rels[1:])): rec(j + 1, cur)
            cur.pop()
    rec(0, []); return out

rows = []
for rl in allLNAs(n):
    alg = nakayama.LinearNakayamaAlgebra(n, rl)
    dual = lm._quiet(pathAlgebra.dualPathAlgebra, alg)
    res = {}
    for v in range(1, n + 1):
        for s in (v, -v):
            if not lm._quiet(mutation.mutationIsPossibleAtVertex, alg if s > 0 else dual, v): continue
            nx_ = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s])
            if nx_ is None: continue
            new = lm.asRelLengths(nx_, n)
            if new is not None: res.setdefault('R' if s > 0 else 'L', []).append((v, new))
    rows.append((rl, res))
    if show: print(rl, res)
cnt = collections.Counter((len(r.get('R', [])), len(r.get('L', []))) for _, r in rows)
print("n", n, "LNAs", len(rows), "(#R,#L) counts", dict(cnt))
print("R v hist", sorted(collections.Counter(r['R'][0][0] for _, r in rows if 'R' in r).items()))
print("L v hist", sorted(collections.Counter(r['L'][0][0] for _, r in rows if 'L' in r).items()))
json.dump(rows, open(f"workshop/rounds/059/experimentalist_which_{n}.json", "w"))
