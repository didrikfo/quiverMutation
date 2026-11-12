"""Round 125 (experimentalist): revision of rounds/123/experimentalist.md after maverick's review.  Repo root:
  timeout 10m .venv/bin/python workshop/rounds/125/experimentalist_transit.py N H [moves ...]     (default: rules mut der)
For the lone-3 key class at (n, h) (rounds/123/experimentalist_classlevel.keyclass; image key = key of m[1:]) and each move family:
 (1) induced moves: every move between two K >= 1 rows (vertex 1 free), is m[1:] -> r[1:] equal or a library move of the
     same family at n - 1?  Control (must say NO somewhere): the same test on moves with one end K >= 1 and the other K = 0.
 (2) restricted components = union-find on K >= 1 rows only; mixed at each K = 1..K0 (maverick_k1.py); shuffled-label control.
 (3) per full orbit: K0 rows by image key; number of restricted components inside it and how many distinct image keys they carry.
 (4) K = 0 transit: components of the K = 0 subgraph; which of them touch restricted components of >= 2 image keys (bridges).
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


def run(n, h, moves, cls, lab):
    t = time.time(); K0 = min(h, n-4-h); K = {m: E.Kof(m, n) for m in cls}
    nb = {m: E.neighbours(n, m, moves) for m in cls}
    edges = [(m, r) for m in cls for r in nb[m]]
    full = uf(cls, edges)
    free = [m for m in cls if K[m] >= 1]; F = set(free)
    # (1) induced moves
    ind = collections.Counter(); ctl = collections.Counter(); imgnb = {}
    def inb(a):
        if a not in imgnb: imgnb[a] = set(E.neighbours(n-1, a, moves))
        return imgnb[a]
    for m, r in edges:
        a, b = m[1:], r[1:]
        bucket = ind if (m in F and r in F) else (ctl if (m in F) != (r in F) else None)
        if bucket is None: continue
        v = 'same image' if a == b else ('image move' if b in inb(a) or a in inb(b) else 'neither')
        bucket[v] += 1
        if bucket is ctl and lab[m] != lab[r]: bucket['(key changes; of these neither: %s)' % (v == 'neither')] += 1
    print(' [%s] induced, both ends K>=1: %s | control, one end K=0: %s' % (moves, dict(ind), dict(ctl)))
    # (2) restricted components
    rc = uf(free, [(m, r) for m, r in edges if m in F and r in F])
    comps = collections.defaultdict(list)
    for m in free: comps[rc[m]].append(m)
    keys_of = {c: {lab[m] for m in v} for c, v in comps.items()}
    line = []
    rnd = random.Random(7)
    for k in range(1, K0+1):
        rowsk = [m for m in free if K[m] == k]
        def nmix(L):
            d = collections.defaultdict(set)
            for m, l in zip(rowsk, L): d[rc[m]].add(l)
            return sum(1 for v in d.values() if len(v) > 1), len(d)
        real, hit = nmix([lab[m] for m in rowsk]); sh = []
        for _ in range(50):
            L = [lab[m] for m in rowsk]; rnd.shuffle(L); sh.append(nmix(L)[0])
        line.append('K=%d hit %d mixed %d (shuffled %.1f, %d-%d)' % (k, hit, real, sum(sh)/50, min(sh), max(sh)))
    print(' [%s] restricted (K>=1) rows %d, components %d, whole component mixed %d; %s; image keys on K>=1 rows: %d'
          % (moves, len(free), len(comps), sum(len(v) > 1 for v in keys_of.values()), '; '.join(line), len({lab[m] for m in free})))
    # (3) per full orbit
    orb = collections.defaultdict(list)
    for m in cls: orb[full[m]].append(m)
    stats = collections.Counter(); rows = []
    for o, v in orb.items():
        k0 = collections.Counter(lab[m] for m in v if K[m] == K0)
        rcs = {rc[m] for m in v if m in F}; ik = set().union(*[keys_of[c] for c in rcs]) if rcs else set()
        if sum(k0.values()) >= 2:
            stats['orbits with >= 2 K0 rows'] += 1
            stats['  of these, > 1 restricted component'] += len(rcs) > 1
            stats['  of these, mixed at K0'] += len(k0) > 1
        stats['orbits joining restricted comps of >= 2 image keys'] += len(ik) > 1
        if len(k0) > 1 or len(ik) > 1 or len(orb) <= 4:
            byk = {k: dict(sorted(collections.Counter(lab[m] for m in v if K[m] == k).items())) for k in range(1, K0+1)}
            rows.append((len(v), len(rcs), len(ik), dict(sorted(k0.items())), byk))
    print(' [%s] full orbits %d; %s' % (moves, len(orb), dict(stats)))
    for r in sorted(rows, reverse=True)[:6]:
        print('    orbit %d rows: %d restricted comps carrying %d image keys; K0 rows by image key %s; K=1..K0 rows by image key %s' % r)
    # (4) K = 0 transit components
    z = [m for m in cls if K[m] == 0]; Z = set(z)
    zc = uf(z, [(m, r) for m, r in edges if m in Z and r in Z])
    touch = collections.defaultdict(set); size = collections.Counter(zc[m] for m in z)
    for m, r in edges:
        for a, b in ((m, r), (r, m)):
            if a in Z and b in F: touch[zc[a]].add(rc[b])
    br = [c for c, t in touch.items() if len(set().union(*[keys_of[x] for x in t])) > 1]
    joins = [c for c, t in touch.items() if len(t) > 1]
    print(' [%s] K=0 rows %d in %d K=0 components; %d touch >= 2 restricted comps; %d touch comps of >= 2 image keys (bridges), sizes %s'
          % (moves, len(z), len(size), len(joins), len(br), sorted((size[c] for c in br), reverse=True)[:8]))
    if br:
        c = min(br, key=lambda c: (size[c], c)); ex = sorted(m for m in z if zc[m] == c)[:3]
        print('    smallest bridge: %d K=0 rows, e.g. %s, touching %d restricted comps' % (size[c], [s(m) for m in ex], len(touch[c])))
    print(' [%s] %.0fs' % (moves, time.time()-t), flush=True)


if __name__ == '__main__':
    n, h = int(sys.argv[1]), int(sys.argv[2]); fams = sys.argv[3:] or ['rules', 'mut', 'der']
    t = time.time(); cls, _ = E.keyclass(n, h); lab = {}; ids = {}
    for m in cls: lab[m] = ids.setdefault(E.key(n-1, m[1:]), len(ids))
    print('n=%d h=%d K0=%d: class %d rows, %d image keys  %.0fs' % (n, h, min(h, n-4-h), len(cls), len(ids), time.time()-t), flush=True)
    for mv in fams: run(n, h, mv, cls, lab)
