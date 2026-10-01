"""At one n: the two orbits P = orbit(word@0), Q = orbit(word@1); which single-word placements lie in them,
and the distribution of row statistics in each.
  timeout 10m .venv/bin/python workshop/rounds/011/theorist_PQ.py N [word]     (word default 35 (even n) / 36 (odd n))
"""
import sys
sys.path.insert(0, '.')
import batch
from collections import Counter
from quivermutation import freeMoves as fm

n = int(sys.argv[1]); w = sys.argv[2] if len(sys.argv) > 2 else ('35' if n % 2 == 0 else '36')
orb = []
for o in (0, 1):
    rep = fm.orbitReport(n, tuple(batch._rowFor(n, w, o)), free=fm.REDUCED, limit=1500000)
    orb.append(rep.rows); print('orbit of', w, '@', o, 'size', len(rep.rows), 'closed', rep.closed, flush=True)
P, Q = orb
print('shared', len(P & Q))
def core(row):
    nz = [i for i, r in enumerate(row) if r]
    if not nz: return None
    # single word = nothing but the word, i.e. whole row is zero-padded word
    return ''.join(map(str, row[nz[0]:nz[-1] + 1])), nz[0]
for name, S in (('P', P), ('Q', Q)):
    words = Counter(); 
    for r in S:
        c = core(r)
        if c and max(r) <= 9: words[c[0]] += 1
    cat = {}
    for r in S:
        c = core(r)
        if c and max(r) <= 9 and len(c[0]) <= 4: cat.setdefault(c[0], []).append(c[1])
    print(name, 'catalogue-like words (len<=4):', {k: sorted(v) for k, v in sorted(cat.items(), key=lambda kv: (len(kv[0]), kv[0]))})
    print(name, 'count dist', sorted(Counter(sum(1 for x in r if x) for r in S).items()))
    print(name, 'sum dist', sorted(Counter(sum(r) for r in S).items()))
    print(name, 'sum mod2 of start positions*len', sorted(Counter(sum((i + 1) * r for i, r in enumerate(r)) % 2 for r in S).items()))
