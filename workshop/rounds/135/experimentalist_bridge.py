"""Round 135 (experimentalist): S-1 / H-026, the K = 0 transit claim of E-325 on lone-3 key classes NOT used in E-325.  Repo root:
  timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_bridge.py N H [moves ...]   (default: rules mut der)
Class = all rows at length n with the Coxeter key of the lone 3 at h (rounds/123/experimentalist_classlevel.keyclass);
K = free head stretch; image key = key(n-1, m[1:]) (head deletion).  Per move family:
 C1  every move edge whose two ends have different image keys: how many have both ends K >= 1 (E-325 says 0),
     one end K = 0 (boundary), both ends K = 0 (inside the K = 0 subgraph).
 C2  per orbit: restricted (K >= 1) components and their image keys; K = 0 components touching restricted comps of >= 2
     keys (bridges), with the largest K reached on each key side of each bridge and the bridge's own internal key changes.
Planted control (the test must say NO): the same C1 count with the label replaced by the TAIL-deletion key key(n-1, m[:-1])
(an invariant for which no induced-move argument exists), and with a random relabelling of the image keys per row.
"""
import sys, os, time, collections, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '123'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '117'))
import experimentalist_classlevel as E


def uf(nodes, edges):
    par = {m: m for m in nodes}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in edges:
        ra, rb = f(a), f(b)
        if ra != rb: par[ra] = rb
    return {m: f(m) for m in nodes}


def s(m): return ''.join(map(str, m))


def c1(edges, K, lab):
    c = collections.Counter()
    for m, r in edges:
        if m == r or lab[m] == lab[r]: continue
        c['both>=1' if K[m] >= 1 and K[r] >= 1 else ('both0' if K[m] == 0 and K[r] == 0 else 'boundary')] += 1
    return dict(c)


def run(n, h, moves, cls, lab, tail, rnd):
    t = time.time(); K0 = min(h, n-4-h); K = {m: E.Kof(m, n) for m in cls}; S = set(cls)
    edges = []
    for m in cls:
        for r in E.neighbours(n, m, moves):
            assert r in S, ('left key class', m, r)
            edges.append((m, r))
    print(' [%s] %d directed edges; C1 key-changing edges by ends: image %s | CONTROL tail-deletion %s | CONTROL random %s'
          % (moves, len(edges), c1(edges, K, lab), c1(edges, K, tail), c1(edges, K, rnd)), flush=True)
    full = uf(cls, edges)
    F = {m for m in cls if K[m] >= 1}; Z = set(cls) - F
    rc = uf(list(F), [(m, r) for m, r in edges if m in F and r in F])
    keys_of = collections.defaultdict(set); maxK = collections.Counter()
    for m in F: keys_of[rc[m]].add(lab[m]); maxK[rc[m]] = max(maxK[rc[m]], K[m])
    zc = uf(list(Z), [(m, r) for m, r in edges if m in Z and r in Z])
    touch = collections.defaultdict(set)
    for m, r in edges:
        for a, b in ((m, r), (r, m)):
            if a in Z and b in F: touch[zc[a]].add(rc[b])
    zsize = collections.Counter(zc[m] for m in Z); zkeys = collections.defaultdict(set)
    for m in Z: zkeys[zc[m]].add(lab[m])
    orb = collections.defaultdict(list)
    for m in cls: orb[full[m]].append(m)
    for o, v in sorted(orb.items(), key=lambda x: -len(x[1])):
        rcs = {rc[m] for m in v if m in F}; ik = set().union(*[keys_of[c] for c in rcs]) if rcs else set()
        k0 = collections.Counter(lab[m] for m in v if K[m] == K0)
        if len(ik) < 2 and len(orb) > 3: continue
        zs = {zc[m] for m in v if m in Z}
        br = [z for z in zs if len(set().union(*[keys_of[c] for c in touch[z]]) if touch[z] else set()) > 1]
        print('   orbit %d rows (K0 rows by image key %s): %d restricted comps carrying %d image keys; %d K=0 comps, %d bridges'
              % (len(v), dict(sorted(k0.items())), len(rcs), len(ik), len(zs), len(br)))
        for z in sorted(br, key=lambda z: -zsize[z])[:5]:
            side = collections.defaultdict(int)
            for c in touch[z]:
                for k in keys_of[c]: side[k] = max(side[k], maxK[c])
            ex = min(m for m in Z if zc[m] == z)
            print('     bridge %d K=0 rows, own image keys %s, touches %d restricted comps; max K per image key %s; e.g. %s'
                  % (zsize[z], sorted(zkeys[z]), len(touch[z]), dict(sorted(side.items())), s(ex)))
    print(' [%s] %.0fs' % (moves, time.time()-t), flush=True)


if __name__ == '__main__':
    n, h = int(sys.argv[1]), int(sys.argv[2]); fams = sys.argv[3:] or ['rules', 'mut', 'der']
    t = time.time(); cls, _ = E.keyclass(n, h); ids = {}; tids = {}; lab = {}; tail = {}
    for m in cls:
        lab[m] = ids.setdefault(E.key(n-1, m[1:]), len(ids))
        tail[m] = tids.setdefault(E.key(n-1, m[:-1]), len(tids))
    rr = random.Random(135); rnd = {m: rr.randrange(len(ids)) for m in cls}
    tab = collections.defaultdict(collections.Counter)
    for m in cls: tab[E.Kof(m, n)][lab[m]] += 1
    print('n=%d h=%d K0=%d: class %d rows, %d image keys, %d tail keys  %.0fs' % (n, h, min(h, n-4-h), len(cls), len(ids), len(tids), time.time()-t))
    print(' image keys by K:', {k: dict(sorted(c.items())) for k, c in sorted(tab.items())}, flush=True)
    for mv in fams: run(n, h, mv, cls, lab, tail, rnd)
