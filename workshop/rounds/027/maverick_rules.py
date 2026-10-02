"""S-1: per same-class set at length n, (a) fixed-position purity, (b) coverage upper bound for any choice rule, (c) structural rules.
usage: python workshop/rounds/027/maverick_rules.py N   (repository root). Class '?' (unresolved cospectral) is dropped, not guessed."""
import sys, os, collections
sys.path.insert(0, os.path.dirname(__file__))
from maverick_classes import classes
from quivermutation import piecewiseHereditary as pwh
n = int(sys.argv[1])
lab, bad = classes(n); lab2, _ = classes(n - 1)
members = collections.defaultdict(list)
for m, l in lab.items():
    if isinstance(l, tuple) and l[1] == '?': continue
    members[l].append(m)
def spans(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]   # relation from vertex s to s+a
def cover(m, v): return sum(1 for s, e in spans(m) if s < v < e)         # relations with v interior
def touch(m, v): return sum(1 for s, e in spans(m) if s <= v <= e)
def image(m, v):
    r = pwh.removeVertex(n, list(m), v)
    return None if r is None else lab2[tuple(r[1])]
RULES = {}
def mk(name, key, pick):
    RULES[name] = (key, pick)
for nm, f in (('cover', cover), ('touch', touch)):
    for pk in ('left', 'right', 'mid'):
        RULES[nm + '-min-' + pk] = (f, pk)
def choose(m, key, pick):
    vs = range(1, n + 1); best = min(key(m, v) for v in vs)
    c = [v for v in vs if key(m, v) == best]
    return c[0] if pick == 'left' else c[-1] if pick == 'right' else min(c, key=lambda v: abs(2 * v - n - 1))
print('n', n, 'classes', len(members), 'LNAs', sum(map(len, members.values())))
print('class key | size | best fixed vertex purity (v) | coverage upper bound | rule purities: ' + ' '.join(RULES))
tot = collections.Counter(); sizes = 0
for c, mem in sorted(members.items(), key=lambda t: -len(t[1])):
    if len(mem) < 2: continue
    pur = {}
    for v in range(1, n + 1):
        cnt = collections.Counter(image(m, v) for m in mem); pur[v] = cnt.most_common(1)[0][1] / len(mem)
    bv = max(pur, key=pur.get)
    sets = [{image(m, v) for v in range(1, n + 1)} - {None} for m in mem]
    allt = collections.Counter(t for s in sets for t in s)
    cov = max(allt.values()) / len(mem)
    rp = []
    for name, (key, pick) in RULES.items():
        cnt = collections.Counter(image(m, choose(m, key, pick)) for m in mem); rp.append(cnt.most_common(1)[0][1] / len(mem))
        tot[name] += cnt.most_common(1)[0][1]
    sizes += len(mem); tot['fixed'] += max(pur.values()) * len(mem); tot['cov'] += cov * len(mem)
    print(str(c)[-34:], len(mem), '%d:%.2f' % (bv, pur[bv]), '%.2f' % cov, ' '.join('%.2f' % x for x in rp))
print('size-weighted over classes of size >= 2:', sizes)
for k, v in tot.items(): print('  %-18s %.3f' % (k, v / sizes))
