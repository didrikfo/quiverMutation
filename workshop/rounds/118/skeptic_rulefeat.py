"""Skeptic r118, item 14: does ANY simple feature of the rule data (old/new window pair lists) predict the delta label (1,1) vs (1,2)?
Pre-registered: an exploratory feature is 'real' only if trained on widths 5+6 (68 rules) it scores >= 90% on width 7 (94 rules) AND beats the
majority-class baseline; features are chosen on 5+6 only.  Control: label-shuffled features must fail the same cut."""
import json, itertools, random, collections
D = {w: json.load(open('workshop/rounds/118/skeptic_ruledelta_w%d.json' % w)) for w in (5, 6, 7)}
def feats(r):
    w, old, new = r
    old = [tuple(x) for x in old]; new = [tuple(x) for x in new]
    f = {'nold': len(old), 'nnew': len(new), 'diff': len(new) - len(old), 'w': w}
    f['maxlen_new'] = max(l for _, l in new); f['maxlen_old'] = max(l for _, l in old)
    f['minlen_new'] = min(l for _, l in new); f['minlen_old'] = min(l for _, l in old)
    f['sumlen_new'] = sum(l for _, l in new); f['sumlen_old'] = sum(l for _, l in old)
    f['span_new'] = max(o + l for o, l in new) - min(o for o, _ in new)
    f['span_old'] = max(o + l for o, l in old) - min(o for o, _ in old)
    f['nnew_minus_nold_pairs_shared'] = len(set(new) - set(old))
    f['maxlen_new_eq_w-1'] = int(f['maxlen_new'] == w - 1); f['maxlen_old_eq_w-1'] = int(f['maxlen_old'] == w - 1)
    f['rel_sumdiff'] = f['sumlen_new'] - f['sumlen_old']
    return f
def rows(ws): return [(feats(r), tuple(y)) for w in ws for r, y in D[w]]
def score(train, test, names):
    tab = collections.defaultdict(collections.Counter)
    for f, y in train: tab[tuple(f[n] for n in names)][y] += 1
    maj = collections.Counter(y for _, y in train).most_common(1)[0][0]
    ok = sum(((tab[tuple(f[n] for n in names)].most_common(1) or [(maj,)])[0][0] if tuple(f[n] for n in names) in tab else maj) == y for f, y in test)
    return ok / len(test)
tr, te = rows((5, 6)), rows((7,))
names = [n for n in tr[0][0] if n != 'w']
base = collections.Counter(y for _, y in te).most_common(1)[0][1] / len(te)
res = []
for k in (1, 2):
    for c in itertools.combinations(names, k): res.append((score(tr, te, c), score(tr, tr, c), c))
res.sort(reverse=True)
print('width-7 majority baseline', round(base, 3), '; features', len(names), '; tried', len(res))
for r in res[:6]: print('  test7 %.3f  train %.3f  %s' % r)
rng = random.Random(0); best = []
for t in range(20):
    ys = [y for _, y in tr]; rng.shuffle(ys); tr2 = [(f, y) for (f, _), y in zip(tr, ys)]
    best.append(max(score(tr2, te, c) for k in (1, 2) for c in itertools.combinations(names, k)))
print('label-shuffled best test7 over', len(best), 'shuffles: max %.3f mean %.3f' % (max(best), sum(best) / len(best)))
