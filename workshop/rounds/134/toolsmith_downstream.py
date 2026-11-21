"""Round 134 (toolsmith), T10, revision of round 133 (E-341): downstream power for the 3 refusals at E-332's roots.

For ROOT (n = 11, F-only library walk, defaults, fingerprint.Visited, depth 6, option OFF) this
  1. records every node of the default ball (canonicalKey) and every gate-admitted key-kept NON-tilting step
     (the steps `refuseNonTilting=True` refuses), with the step number at which it occurs;
  2. builds each refused child (as the library builds it) and checks it is in the default ball;
  3. from each refused child runs an F-only depth-K walk (option off, key guard on with the ROOT's key) and reports
     the nodes it reaches that are NOT in the root's default depth-6 ball, split in-class (same SNF(C + C^T) and
     signature as the root) / cross-class; an in-class such node is a candidate node the option would lose at depth 6 + K;
  4. CONTROL (detector power): the same depth-K walk from up to M in-class TILTING-step children at step 6
     (nodes the option keeps); it must report in-class nodes outside the ball (else the detector cannot say YES).
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/134/toolsmith_downstream.py 303500030 [root|dual] [K=3] [M=2]
"""
import sys, os, json, time, copy, collections, importlib.util
sys.path.insert(0, '.')
spec = importlib.util.spec_from_file_location('dcount', 'workshop/rounds/133/toolsmith_dcount.py')
dc = importlib.util.module_from_spec(spec); spec.loader.exec_module(dc)   # wraps the gate with its counter (inactive unless ST['on'])
from quivermutation import search, mutation, fingerprint

def ball(start, D, base, refuse=False):
    got = {}
    def visit(A, path): got.setdefault(fingerprint.canonicalKey(A), (copy.deepcopy(A), list(path)))
    search.mutationSearchDepthFirst(copy.deepcopy(start), D, printOutput=False, visitor=visit, visited=fingerprint.Visited(),
                                    baseKey=base, refuseNonTilting=refuse)
    return got

def main(root, side, K, M):
    t0 = time.time(); N, D = 11, 6
    s = dc.startOf(N, root, side); i0 = dc.inv(s); base = search._coxeterKeyOrNone(s)
    # 1. default ball, with the step depth of every gate call (depth = len(path of parent) + 1)
    steps = []; depthOf = {}
    gate = mutation.mutationIsPossibleAtVertex
    got = {}
    def visit(A, path):
        k = fingerprint.canonicalKey(A); got.setdefault(k, (copy.deepcopy(A), list(path))); depthOf[id(A)] = len(path)
    def spy(alg, v):
        ok = gate(alg, v)
        if ok and id(alg) in depthOf:
            t = dc._tilt(alg, v)
            if (not t) or (depthOf[id(alg)] + 1 == D and sum(r['tilting'] for r in steps) < 10 * M):
                legal, c = dc.childOf(alg, v)
                if legal and search._coxeterKeyOrNone(c) == base:
                    steps.append(dict(step=depthOf[id(alg)] + 1, v=v, tilting=t, child=c))
        return ok
    mutation.mutationIsPossibleAtVertex = spy
    search.mutationSearchDepthFirst(copy.deepcopy(s), D, printOutput=False, visitor=visit, visited=fingerprint.Visited(), baseKey=base)
    mutation.mutationIsPossibleAtVertex = gate
    ballKeys = set(got); tBall = round(time.time() - t0)
    refused = [r for r in steps if not r['tilting']]
    seenC = set(); uniq = []
    for r in refused:
        k = fingerprint.canonicalKey(r['child'])
        if k not in seenC: seenC.add(k); uniq.append(r)
    print(json.dumps(dict(root=root, side=side, ballNodes=len(ballKeys), secsBall=tBall, refusedKeptSteps=len(refused),
                          refusedAtStep=collections.Counter(r['step'] for r in refused), distinctRefusedChildren=len(uniq))), flush=True)
    # 3. downstream walks
    def downstream(r, tag):
        c = r['child']; ck = fingerprint.canonicalKey(c); t1 = time.time()
        dc.ST.update(on=True, base=base, C=collections.Counter(), rows=[])
        g = ball(c, K, base); dc.ST['on'] = False
        out = [k for k in g if k not in ballKeys]
        lab = collections.Counter('in' if dc.inv(g[k][0]) == i0 else 'cross' for k in out)
        inPaths = [' '.join(map(str, g[k][1])) for k in out if dc.inv(g[k][0]) == i0][:5]
        print(json.dumps(dict(tag=tag, step=r['step'], v=r['v'], childInBall=ck in ballKeys, childClass='in' if dc.inv(c) == i0 else 'cross',
                              childSnfTail=dc.inv(c)[0][-3:], nodes=len(g), outsideBall=len(out), outsideByClass=dict(lab),
                              inClassOutsidePaths=inPaths, downstreamNonTiltingKept=len(dc.ST['rows']),
                              downstreamCounts={k: v for k, v in dc.ST['C'].items() if k in ('admitted', 'd>=2', 'notTilting', 'NONVACUOUS')},
                              secs=round(time.time() - t1))), flush=True)
        return lab.get('in', 0)
    lossCand = sum(downstream(r, 'REFUSED') for r in uniq)
    ctl = []; seenT = set()
    for r in steps:
        if not r['tilting'] or r['step'] != D: continue
        k = fingerprint.canonicalKey(r['child'])
        if k in seenT or dc.inv(r['child']) != i0: continue
        seenT.add(k); ctl.append(r)
        if len(ctl) >= M: break
    ctlHits = [downstream(r, 'CONTROL tilting step-6 child') for r in ctl]
    print('SUMMARY root %s %s: refused children %d, in-class nodes outside the depth-6 ball reached from them at +%d: %d; '
          'control children %d, in-class outside-ball nodes: %s (control %s); %.0fs'
          % (root, side, len(uniq), K, lossCand, len(ctl), ctlHits, 'PASS' if ctl and all(h > 0 for h in ctlHits) else 'FAIL', time.time() - t0), flush=True)

if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1] if len(a) > 1 else 'root', int(a[2]) if len(a) > 2 else 3, int(a[3]) if len(a) > 3 else 2)
