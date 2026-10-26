"""Round 089 theorist, T10/T5: can Oppermann's mutation (arXiv:1504.02617v1, Sect. 3, Steps 1-4) at a
vertex i produce a loop of degree >= 1 at i?  Quiver level only (arrows and degrees, no differential).
Arrows: (name, s, t, deg).  Step 1 rotation, Step 2 alpha.phi, Step 3 phi.alpha^-1, Step 4 alpha.phi.beta^-1.
Part A: one step from random ungraded parents (parallel arrows, cycles allowed) -> which M-arrows are loops, and
        compare with the predicted list (a) loop of Q (b) 2-cycle through i (c) parallel i->v with a deg>=1 arrow
        (d) deg>=1 loop at i (Step 4).
Part B: iterated mutation sequences (depth <= D) from random ACYCLIC ungraded parents, n <= 5: does any loop
        (any vertex, any degree) ever occur, and is the hypothesis 'no degree-0 loop at i' ever violated?
Usage: .venv/bin/python workshop/rounds/089/theorist_loopreach.py [ntrials] [depth]
"""
import sys, random, itertools

def mutate(arrows, i):
    A = [a for a in arrows if a[1] == i and a[3] == 0]
    assert all(a[2] != i for a in A), "degree-0 loop at i"
    Aname = {a[0] for a in A}
    out = []
    for (nm, s, t, d) in arrows:
        if nm in Aname: continue
        out.append((nm, s, t, d - (s == i) + (t == i)))                 # step 1
    for (nm, s, t, d) in A: out.append((('*', nm), t, i, 0))
    for (nm, s, t, d) in arrows:                                         # step 2
        if t == i and nm not in Aname:
            for a in A: out.append((('o', a[0], nm), s, a[2], d - (s == i) + 1 - 1))
    for (nm, s, t, d) in arrows:                                         # step 3
        if s == i and nm not in Aname:
            for a in A: out.append((('i', nm, a[0]), a[2], t, d - 1 + 1 + (t == i) * 0 + 0 if False else (d - 1 + (t == i)) + 1))
    for (nm, s, t, d) in arrows:                                         # step 4
        if s == i and t == i and nm not in Aname:
            for a in A:
                for b in A: out.append((('c', a[0], nm, b[0]), b[2], a[2], d))
    return out

def loops(arrows): return [a for a in arrows if a[1] == a[2]]

def rand_q(rng, n, acyclic, graded):
    ar = []; k = 0
    for _ in range(rng.randint(1, 3 * n)):
        s, t = rng.randrange(n), rng.randrange(n)
        if acyclic:
            if s == t: continue
            s, t = min(s, t), max(s, t)
        if s == t and not graded: continue
        ar.append((('q', k), s, t, rng.choice([1, 2]) if (graded and rng.random() < .4) else 0)); k += 1
    return ar

def has_cycle(arrows, n):
    adj = {v: set() for v in range(n)}
    for _, s, t, _ in arrows: adj[s].add(t)
    seen = {}
    def dfs(v):
        seen[v] = 1
        for w in adj[v]:
            if seen.get(w) == 1 or (w not in seen and dfs(w)): return True
        seen[v] = 2; return False
    return any(v not in seen and dfs(v) for v in range(n))

def partA(ntr, seed=1):
    rng = random.Random(seed); n_steps = 0; bad = 0; kinds = {}; deg_ge1_at_i = 0
    for _ in range(ntr):
        n = rng.randint(2, 5); Q = rand_q(rng, n, False, rng.random() < .5)
        i = rng.randrange(n)
        if any(a[1] == i == a[2] and a[3] == 0 for a in Q): continue
        M = mutate(Q, i); n_steps += 1
        A = {a[0] for a in Q if a[1] == i and a[3] == 0}
        for L in loops(M):
            nm = L[0]; v = L[1]; why = None
            if not isinstance(nm, tuple) or nm[0] == 'q': why = 'a:loop of Q'
            elif nm[0] == 'o': why = 'b:2-cycle'
            elif nm[0] == 'i': why = 'c:parallel i->v with deg>=1' if True else None
            elif nm[0] == 'c': why = 'd:step4'
            kinds[why] = kinds.get(why, 0) + 1
            if v == i and L[3] >= 1: deg_ge1_at_i += 1
            # prediction check
            if why[0] == 'a': ok = any(q[0] == nm for q in Q)
            elif why[0] == 'b': ok = Q_has_2cycle(Q, nm, i, v)
            elif why[0] == 'c': ok = (nm[2] in A) and any(q[0] == nm[1] and q[3] >= 1 and q[1] == i for q in Q)
            else: ok = any(q[0] == nm[2] and q[1] == q[2] == i and q[3] >= 1 for q in Q)
            if not ok: bad += 1
    return n_steps, kinds, bad, deg_ge1_at_i

def Q_has_2cycle(Q, nm, i, v):
    return any(q[0] == nm[2] and q[1] == v and q[2] == i for q in Q) and any(q[0] == nm[1] and q[1] == i and q[2] == v for q in Q)

def partB(ntr, depth, seed=2):
    rng = random.Random(seed); seqs = 0; first_loop = 0; viol = 0; maxdeg = 0; stat = {}
    for _ in range(ntr):
        n = rng.randint(2, 5); Q = rand_q(rng, n, True, False)
        if not Q: continue
        def rec(Q, d):
            nonlocal seqs, first_loop, viol, maxdeg
            seqs += 1
            maxdeg = max([maxdeg] + [a[3] for a in Q])
            if loops(Q):
                first_loop += 1
                for L in loops(Q):
                    key = (d, 'deg%d' % L[3]); stat[key] = stat.get(key, 0) + 1
            if d == depth: return
            for i in range(n):
                if any(a[1] == i == a[2] and a[3] == 0 for a in Q): viol += 1; continue
                rec(mutate(Q, i), d + 1)
        rec(Q, 0)
    print('  loops by (depth, degree):', sorted(stat.items()))
    return seqs, first_loop, viol, maxdeg

if __name__ == '__main__':
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    depth = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    print('A one step, random parents (cycles, parallel arrows, degrees 0..2 allowed):', partA(ntr))
    print('B iterated, acyclic ungraded start: (nodes, nodes with a loop, steps blocked by deg-0 loop at i, max arrow degree):', partB(max(ntr // 10, 50), depth))

def minimal_example(maxn=4, norepeat=False):
    """smallest acyclic ungraded parent (fewest vertices, then arrows) and a shortest mutation sequence ending in a loop"""
    best = None
    for n in range(2, maxn + 1):
        pairs = [(s, t) for s in range(n) for t in range(s + 1, n)]
        for m in range(1, 5):
            for combo in itertools.combinations_with_replacement(pairs, m):
                Q = [(('q', k), s, t, 0) for k, (s, t) in enumerate(combo)]
                stack = [(Q, [])]
                for depth in range(5):
                    nxt = []
                    for (R, seq) in stack:
                        for i in range(n):
                            if norepeat and seq and seq[-1] == i: continue
                            if any(a[1] == i == a[2] and a[3] == 0 for a in R): continue
                            R2 = mutate(R, i)
                            if loops(R2): return n, combo, seq + [i], loops(R2)[0], i
                            nxt.append((R2, seq + [i]))
                    stack = nxt
    return best
