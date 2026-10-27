"""Round 090 theorist, T10/T5 (revision of round 089).  Usage (repo root):
  .venv/bin/python workshop/rounds/090/theorist_relarrows.py ex                 # the omega counterexample
  .venv/bin/python workshop/rounds/090/theorist_relarrows.py rel [ntrials] [depth]   # Parts A/B with relation arrows
  .venv/bin/python workshop/rounds/090/theorist_relarrows.py parents /tmp/tsm/c1.pkl 1   # 25 failing parents: acyclic with relations
Parents for 'rel': random ACYCLIC quiver (degree-0 arrows, parallel allowed) plus, for random pairs s<t joined by a
degree-0 path of length >= 2, a relation arrow omega: s->t of degree 1 (Oppermann const.usual_algebra, l. 291), and
sometimes a degree-2 arrow s->t (syzygy of the omega; Ext^2 between s and t needs a path s->t too).
Quiver level only (arrows, endpoints, degrees)."""
import sys, random, itertools
sys.path.insert(0, 'workshop/rounds/089')
from theorist_loopreach import mutate, loops

def reach_len2(arrows, s, t):
    """is there a degree-0 path s -> t of length >= 2 ?"""
    out = {}
    for a in arrows:
        if a[3] == 0: out.setdefault(a[1], []).append(a[2])
    return any(t in out.get(m, []) or any(t in out.get(k, []) for k in out.get(m, [])) or _r(out, m, t) for m in out.get(s, []))
def _r(out, m, t):
    seen, st = set(), [m]
    while st:
        x = st.pop()
        for y in out.get(x, []):
            if y == t: return True
            if y not in seen: seen.add(y); st.append(y)
    return False

def rand_parent(rng, n):
    ar = []; k = 0
    for _ in range(rng.randint(1, 2 * n)):
        s, t = sorted(rng.sample(range(n), 2)); ar.append((('q', k), s, t, 0)); k += 1
    rel = []
    for s in range(n):
        for t in range(s + 1, n):
            if reach_len2(ar, s, t) and rng.random() < .5:
                rel.append((('w', k), s, t, 1)); k += 1
                if rng.random() < .3: rel.append((('v', k), s, t, 2)); k += 1
    return ar + rel

def predicted_c_loops(Q, i):
    A = [a for a in Q if a[1] == i and a[3] == 0]
    return sorted((('i', p[0], a[0]), a[2]) for p in Q if p[1] == i and p[3] >= 1 for a in A if a[2] == p[2])

def legal(Q, i): return not any(a[1] == i == a[2] and a[3] == 0 for a in Q)

def main_rel(ntr, depth):
    rng = random.Random(7); steps = bad = at_i = with_loop = nloops = 0; degs = {}
    for _ in range(ntr):
        n = rng.randint(3, 6); Q = rand_parent(rng, n)
        for i in range(n):
            M = mutate(Q, i); steps += 1
            L = loops(M)
            got = sorted((a[0], a[1]) for a in L)
            if got != predicted_c_loops(Q, i): bad += 1
            at_i += any(a[1] == i for a in L); with_loop += bool(L); nloops += len(L)
            for a in L: degs[a[3]] = degs.get(a[3], 0) + 1
    print('A one step, %d parents -> %d steps: loops at i: %d; steps with a loop anywhere: %d (%d loops, degrees %s); '
          'loop set != prediction "pairs (omega out of i, alpha in A, same target)": %d' % (ntr, steps, at_i, with_loop, nloops, degs, bad))
    # Part B: iterate, never repeating the vertex just mutated; count loops at i of the NEXT step's vertex, blocked steps
    rng = random.Random(8); nodes = lp = blocked = posloop = 0; bydepth = {}
    def rec(Q, n, d, last):
        nonlocal nodes, lp, blocked, posloop
        nodes += 1
        if loops(Q): lp += 1; bydepth[d] = bydepth.get(d, 0) + 1
        if d == depth: return
        for i in range(n):
            if i == last: continue
            if not legal(Q, i): blocked += 1; continue
            posloop += any(a[1] == i == a[2] for a in Q); rec(mutate(Q, i), n, d + 1, i)
    for _ in range(max(ntr // 20, 20)):
        n = rng.randint(3, 5); rec(rand_parent(rng, n), n, 0, -1)
    print('B iterated depth <= %d: nodes %d, nodes with a loop %d (by depth %s), steps blocked by a degree-0 loop at i: %d, legal steps at a vertex with a loop of degree >= 1 (E-235 open case): %d' % (depth, nodes, lp, sorted(bydepth.items()), blocked, posloop))

def example():
    Q = [('a', 0, 1, 0), ('b', 1, 2, 0), ('c', 0, 2, 0), ('w', 0, 2, 1)]
    print('referee counterexample (0->1->2, 0->2, relation 0->1->2, omega 0->2 deg 1), mutate at 0:', loops(mutate(Q, 0)))
    Q2 = [('a', 0, 1, 0), ('b', 1, 2, 0), ('w', 0, 2, 1)]
    print('same without the arrow c (A_3 / rad^2 alone, Oppermann l. 296-307), mutate at 0:', loops(mutate(Q2, 0)), ' at 1:', loops(mutate(Q2, 1)), ' at 2:', loops(mutate(Q2, 2)))

def parents(pk, cls):
    import pickle
    from quivermutation import arrowPaths as ap
    recs = [r for r in pickle.load(open(pk, 'rb')) if r['kind'] == 'fail']
    cyc = loopv = relloop = norel = 0; nrel = 0; maxlen = []
    for r in recs:
        alg = r['parentObj']; q = alg.quiver; v = r['v']
        cyc += ap.hasOrientedCycle(q)                              # oriented cycle or loop in the Gabriel quiver
        loopv += any(a == b for a, b in q.edges())
        for rel in alg.rels:                                       # relation arrow omega: first vertex -> last vertex
            nrel += 1; relloop += all(p[0] == p[-1] for p in rel) or rel[0][0] == rel[0][-1]
        norel += not alg.rels
    print('class', cls, 'failing parents', len(recs), '| Gabriel quiver with oriented cycle or loop:', cyc, '| loop edges:', loopv,
          '| relations total', nrel, 'with equal endpoints (omega a loop):', relloop, '| parents without relations:', norel)

if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'ex': example()
    elif m == 'rel': main_rel(int(sys.argv[2]) if len(sys.argv) > 2 else 2000, int(sys.argv[3]) if len(sys.argv) > 3 else 3)
    else: parents(sys.argv[2], int(sys.argv[3]))
