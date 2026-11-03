"""Round 111 theorist: E-240 Lemma A check re-run with arrow degrees in -2..2 (negative degrees included), as the referee asked.
Same engine as workshop/rounds/086/theorist_lemmaA.py; only rand_instance is replaced. At least one arrow of degree -1 is forced.
Usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/111/theorist_lemmaA_neg.py [trials_per_n] [WMAX] [seed]"""
import importlib.util, sys, random
sp = importlib.util.spec_from_file_location('la', 'workshop/rounds/086/theorist_lemmaA.py')
la = importlib.util.module_from_spec(sp); sp.loader.exec_module(la)
def rand_instance(rng, n):
    na = rng.randint(2, 6)
    arrows = [(rng.randrange(n), rng.randrange(n), rng.choice([-2, -1, -1, 0, 0, 1, 2])) for _ in range(na)]
    k = rng.randrange(na); arrows[k] = (arrows[k][0], arrows[k][1], -1)
    idx = list(range(na)); rng.shuffle(idx)
    c = rng.randint(1, na - 1)
    X1 = set(idx[:c]); Y = set(idx[c:]) if rng.random() < .7 else set(idx[c:c + 1])
    if rng.random() < .5 and len(Y) > 1: Y.pop()
    return arrows, sorted(X1), sorted(Y)
la.rand_instance = rand_instance
_rk = la.rank_mod
# (-1)**negative is a float in Python: the sign (-1)**pre must be an exact integer, so cast (values are +-1)
la.rank_mod = lambda rows: _rk([{c: int(round(v)) for c, v in r.items()} for r in rows])
T = int(sys.argv[1]) if len(sys.argv) > 1 else 40
W = int(sys.argv[2]) if len(sys.argv) > 2 else 6
rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
tot = [0, 0, 0, 0, 0]
for n in range(1, 6):
    c = b = tb = neg = 0
    for _ in range(T):
        arrows, X1, Y = rand_instance(rng, n)
        cc, bb, t = la.check(arrows, X1, Y, W)
        c += cc; b += bb; tb += t; neg += sum(1 for a in arrows if a[2] < 0)
    print(f'n={n}: instances={T} negative-degree arrows={neg} cells={c} mismatches={b} overlapping-sites={tb}')
    tot = [tot[0]+T, tot[1]+c, tot[2]+b, tot[3]+tb, tot[4]+neg]
print('TOTAL instances, cells, mismatches, overlaps, neg arrows:', tot)
cb = ci = 0
for _ in range(4 * T):
    arrows, X1, Y = rand_instance(rng, rng.randint(1, 3))
    z = rng.choice(X1); Y = sorted(set(Y) | {z})
    if not any(arrows[x][1] == arrows[y][0] for x in X1 for y in Y): continue
    try: cc, bb, t = la.check(arrows, X1, Y, W)
    except AssertionError: bb = 1
    ci += 1; cb += (bb > 0)
print(f'CONTROL (shared arrow): instances={ci}, mismatching={cb}')
