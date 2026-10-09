"""Round 038 (toolsmith): meet-in-the-middle for E-134's 16 key coincidences. Key-preserving BFS from the 2 LNAs of class 0 and (separately) from
the 16 hit fans, SECONDS for the LNA side, SECONDS/8 per hit; report the number of canonical keys in the intersection (nonzero would put a hit in the LNA class), and the (d,J) of the hit side.
Usage: toolsmith_n6meet.py SECONDS"""
import sys, time, json
secs = float(sys.argv[1]); sys.argv = ['x']
src = open('workshop/rounds/038/toolsmith_n6close.py').read().split("budget = args")[0]
sys.argv = ['x', '--plan', '1']; exec(compile(src, 'c', 'exec'))
key0 = order[0]
def bfs(starts, secs):
    seen = {}; fr = []; t0 = time.time(); tab = Counter()
    def add(x):
        ck = fingerprint.canonicalKey(x)
        if ck in seen: return
        seen[ck] = 1; fr.append(x)
    for x in starts: add(x)
    done = False
    while fr and not done:
        cur, fr[:] = fr[:], []
        for i, alg in enumerate(cur):
            if time.time() - t0 > secs: done = True; break
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                for dj in dims(alg, v).values(): tab[dj] += 1
                raw = mutation.quiverMutationAtVertex(alg, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                ch = reduction.reducePathAlgebra(raw)
                if search._coxeterKeyOrNone(ch) == key0: add(ch)
    return seen, tab, not fr
starts = [build([tuple(x) for x in h[0]], h[1], list(range(1, 7))) for h in hits]
print('hit keys == class-0 key:', all(search._coxeterKeyOrNone(x) == key0 for x in starts), '| hit canonical keys distinct:', len(hitkeys))
S1, t1, c1 = bfs(classes[key0], secs); print('LNA side seen', len(S1), 'closed', c1, '(d,J)', dict(sorted(t1.items())), flush=True)
print('hits inside the LNA-side forward set:', sum(1 for k in hitkeys if k in S1))
tot = set()
for j, st in enumerate(starts):
    S2, t2, c2 = bfs([st], secs / 8); m = set(S1) & set(S2); tot |= m
    print('hit %d (v=%d): seen %d closed %s, meets LNA side in %d algebras | (d,J) with d>=3: %s' % (j, hits[j][3], len(S2), c2, len(m), {k: v for k, v in sorted(t2.items()) if k[0] >= 3}), flush=True)
print('INTERSECTION (union over hits)', len(tot))
