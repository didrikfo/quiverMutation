"""Referee check: omega = ||w|| + stars is non-increasing under d, strictly except the star term; ||.||>=0; multiple seeds/sizes.
Usage: .venv/bin/python workshop/rounds/083/experimentalist_omega.py NTRIALS"""
import sys, random
sys.path.insert(0, '.')
import importlib.util, os
HERE = os.path.dirname(os.path.abspath(__file__))
def load(p, n):
    s = importlib.util.spec_from_file_location(n, os.path.join(HERE, p)); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
sf = load('theorist_starfilt.py', 'sf'); opp = sf.opp
def om(m, w): return sum(m.deg((x,)[0]) for x in w) + sf.nstar(w)
tot = dict(letters=0, neg=0, up=0, eq_nonstar=0, multiA=0)
for seed in range(int(sys.argv[1])):
    rng = random.Random(1000 + seed)
    q = sf.random_parent(rng, False)
    if len(q.A) > 1: tot['multiA'] += 1
    m = opp.M(q, {})
    for L in opp.all_letters(q):
        tot['letters'] += 1
        w0 = (L,)
        if om(m, w0) - sf.nstar(w0) < 0: tot['neg'] += 1
        dl = m.partial({w0: 1})
        for w in dl:
            if om(m, w) > om(m, w0): tot['up'] += 1
            if om(m, w) == om(m, w0) and sf.nstar(w) == 0: tot['eq_nonstar'] += 1
print(tot)
