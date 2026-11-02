"""Round 109 response: per step on the shortest paths, compare repo reduced result with End(T') by canonicalKey (labelled, relation-sensitive), not only Cartan.  Usage: n depth"""
import sys
src = open('workshop/rounds/109/toolsmith_pathaudit.py').read()
exec(src.split("\nif __name__ == '__main__':")[0])
from quivermutation import fingerprint
def ck(a): return fingerprint.canonicalKey(a)
tot = 0; cartan = 0; eq = {'id': 0, 'dual': 0}; bad = 0
for k in hits:
    for pth in paths(k):
        for step in range(len(pth) - 1):
            nxtk = pth[step + 1]; slot = [i for pk, i in parents[nxtk] if pk == pth[step]][0]
            Tn = ordered[pth[step]]; newobj = [x for x in ordered[nxtk] if x not in Tn][0]; T2 = Tn[:slot] + (newobj,) + Tn[slot + 1:]
            pres, _ = presentation(Tn); pres2, _ = presentation(T2)
            a = repoAlg(pres, False); target = repoAlg(pres2, False)
            H2 = Hm(T2); ok = None
            for side in ('right', 'left'):
                b = pathAlgebra.dualPathAlgebra(a) if side == 'left' else a
                if mutation.mutationIsPossibleAtVertex(b, slot + 1) is not True: continue
                f = mutation.leftQuiverMutationAtVertex if side == 'left' else mutation.quiverMutationAtVertex
                red = quiet(reduction.reducePathAlgebra, quiet(f, a, slot + 1)); C = cartanOf(red)
                if not ((C == H2).all() or (C == H2.T).all()): continue
                tot += 1
                target = quiet(reduction.reducePathAlgebra, target); kt = ck(target); kd = ck(pathAlgebra.dualPathAlgebra(target))
                kr = ck(red)
                if kr == kt: eq['id'] += 1
                elif kr == kd: eq['dual'] += 1
                else:
                    bad += 1; print('RELATION MISMATCH step', step, 'slot', slot, side)
                    if bad <= 2: print(' repo  ', kr); print(' target', kt)
print('n', n, 'Cartan-matched steps', tot, 'canonicalKey equal to End(T\'):', eq, 'unequal', bad)
