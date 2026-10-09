"""L1 test: for every LNA of length n (n in range given), every vertex v and direction
(+v right, -v left) with mutationIsPossibleAtVertex: if the result is an LNA, record the
positions j (start vertex of a relation, 1-based) at which the relation-length vector
changed.  Reports, per n, the number of mutation events, how many give LNAs, and the
histogram of offsets (j - v) over changed starts, plus of (end change).  Usage: n_lo n_hi"""
import sys, collections, time
from quivermutation import lnaMoves as lm, mutation, nakayama, pathAlgebra
n_lo, n_hi = int(sys.argv[1]), int(sys.argv[2])

def allLNAs(n):
    out = []
    def rec(j, cur):
        if j == n - 2:
            if lm.isAdmissible(n, cur): out.append(list(cur))
            return
        for a in [0] + list(range(2, n - 1 - j + 1)):
            cur.append(a)
            # prune: starts/ends strictly increase among nonzero
            rels = lm.relationsOf(cur)
            ok = all(e[0] + e[1] < l[0] + l[1] for e, l in zip(rels, rels[1:]))
            if ok: rec(j + 1, cur)
            cur.pop()
    rec(0, [])
    return out

for n in range(n_lo, n_hi + 1):
    t0 = time.time()
    lnas = allLNAs(n)
    ev = lna_ct = 0
    startOff = collections.Counter(); anyChange = collections.Counter(); worst = []
    for rl in lnas:
        alg = nakayama.LinearNakayamaAlgebra(n, rl)
        dual = lm._quiet(pathAlgebra.dualPathAlgebra, alg)
        for v in range(1, n + 1):
            for signed in (v, -v):
                tgt = alg if signed > 0 else dual
                if not lm._quiet(mutation.mutationIsPossibleAtVertex, tgt, v): continue
                ev += 1
                nxt = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [signed])
                if nxt is None: continue
                new = lm.asRelLengths(nxt, n)
                if new is None: continue
                lna_ct += 1
                ch = [j + 1 for j in range(n - 2) if rl[j] != new[j]]
                # a relation change at start j means relation starting at j changed; measure vs v
                for j in ch:
                    startOff[j - v] += 1
                if ch:
                    d = max(abs(j - v) for j in ch)
                    anyChange[d] += 1
                    if d > 2 and len(worst) < 5: worst.append((rl, signed, new))
    print(f"n={n} LNAs={len(lnas)} events={ev} resultIsLNA={lna_ct} ({time.time()-t0:.0f}s)")
    print("  max |start - v| over changed starts:", dict(sorted(anyChange.items())))
    print("  changed-start offset (start - v):", dict(sorted(startOff.items())))
    for w in worst: print("  far example", w)
