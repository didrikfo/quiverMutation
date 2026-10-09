"""Overhead of checkCartan on mutateAtVertex (round 022). Times the 10 E-084 parents' step and a hot-path sample:
mutateAtVertex alone vs with checkCartan=True, then a depth-5 guarded n = 8 sample of admitted steps."""
import sys, time
sys.argv = ['x']
exec(open('workshop/rounds/022/toolsmith_replay.py').read().split("kept = 0")[0])
def time_step(alg, v, k, reps=5):
    rl = procedure.relationsFrom(alg); best = 1e9
    for _ in range(reps):
        t = time.perf_counter(); procedure.mutateAtVertex(alg.quiver, rl, v, checkCartan=k); best = min(best, time.perf_counter() - t)
    return best
tot0 = tot1 = 0
for d, rels, v, path in lines:
    alg = classes[base][path[0]]
    for w in path[1:]:
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, w))
    a, b = time_step(alg, v, False), time_step(alg, v, True)
    tot0 += a; tot1 += b
    print('v', v, 'plain ms %.1f' % (1000 * a), 'checked ms %.1f' % (1000 * b), 'ratio %.2f' % (b / a))
print('10 parents: plain %.1f ms, checked %.1f ms, ratio %.2f' % (1000 * tot0, 1000 * tot1, tot1 / tot0))
# whole step as the walker does it: mutate + reduce + canonicalKey, per expansion
# context: the rest of a walker step (gate + quiverMutationAtVertex + reduce + key) on the same parents
tot2 = 0
for d, rels, v, path in lines:
    alg = classes[base][path[0]]
    for w in path[1:]:
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, w))
    t = time.perf_counter()
    mutation.mutationIsPossibleAtVertex(alg, v)
    ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, v)); search._coxeterKeyOrNone(ch)
    tot2 += time.perf_counter() - t
print('walker-style step (gate + mutate + reduce + key) on the 10 parents: %.1f ms; check adds %.1f ms (+%.0f%%)' % (1000 * tot2, 1000 * (tot1 - tot0), 100 * (tot1 - tot0) / tot2))
