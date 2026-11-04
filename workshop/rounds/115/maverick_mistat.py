"""T6/H-017 revision: mutual information key vs defect (n=9 depth 5, 645 members), several nulls, LNA-cluster rank, planted controls.
usage: .venv/bin/python workshop/rounds/115/maverick_mistat.py   (reads workshop/rounds/105/experimentalist_defect_*.json)"""
import json, glob, collections, random, math
random.seed(7)
rows = []
for f in sorted(glob.glob('workshop/rounds/105/experimentalist_defect_*.json')):
    d = json.load(open(f))
    for m in d['members']:
        rows.append(dict(lna=d['lna'], key=tuple(m['key']), cords=m['cords'], len=m['len'], defect=m['rels'] - m['cords']))
def mi(x, y):
    n = len(x); cx = collections.Counter(x); cy = collections.Counter(y); cxy = collections.Counter(zip(x, y))
    return sum(c / n * math.log(c * n / (cx[a] * cy[b])) for (a, b), c in cxy.items())
def shuf_within(y, cells):
    idx = collections.defaultdict(list)
    for i, c in enumerate(cells): idx[c].append(i)
    y2 = list(y)
    for ii in idx.values():
        v = [y[i] for i in ii]; random.shuffle(v)
        for i, t in zip(ii, v): y2[i] = t
    return y2
def pval(x, y, cells, B=2000):
    o = mi(x, y); ge = 0
    for _ in range(B): ge += mi(x, shuf_within(y, cells)) >= o
    return o, (ge + 1) / (B + 1)
def report(rs, label, y=None):
    y = y or [r['defect'] for r in rs]; x = [r['key'] for r in rs]
    print(label)
    for nm, cf in [('global shuffle', lambda r: 0), ('within (cords,len)', lambda r: (r['cords'], r['len'])), ('within cords', lambda r: r['cords'])]:
        o, p = pval(x, y, [cf(r) for r in rs]); print("  MI=%.4f nats  null: %-20s p=%.4f" % (o, nm, p))
report(rows, "OBSERVED (645)")
# LNA-cluster calibration: key class 2 == the single LNA 3033030. MI of 'LNA a vs rest' with defect, for each LNA a.
print("one-vs-rest MI(defect) per LNA, plus within-(cords,len) p")
y = [r['defect'] for r in rows]; cells = [(r['cords'], r['len']) for r in rows]
for a in sorted(set(r['lna'] for r in rows)):
    x = [r['lna'] == a for r in rows]; o, p = pval(x, y, cells, 500)
    print("  %s n=%3d MI=%.4f p=%.3f" % (a, sum(x), o, p))
# planted controls (same code path as observed)
cls = lambda r: r['key'][2] != 0
d1 = [r['defect'] for r in rows if not cls(r)]; d2 = [r['defect'] for r in rows if cls(r)]
pool = [r['defect'] for r in rows]
def planted(draw, label): 
    yy = [draw(r) for r in rows]; report(rows, label, yy)
planted(lambda r: random.choice(d2) if cls(r) else random.choice(d1), "PLANTED mode-preserving shift (class draws from observed class distributions; YES expected for MI)")
planted(lambda r: random.choice(pool), "PLANTED null: iid from pooled distribution (NO expected)")
planted(lambda r: random.choice([2,3,4]) if cls(r) else random.choice([1,2,2,2,3,3,4,4]), "PLANTED mild shift, same mode 2 (YES expected)")
# purity on the mode-preserving plant: blindness check
def purity(cells, y):
    c = collections.defaultdict(collections.Counter)
    for k, v in zip(cells, y): c[k][v] += 1
    return sum(max(x.values()) for x in c.values()) / len(y)
yy = [random.choice(d2) if cls(r) else random.choice(d1) for r in rows]
kc = [r['key'] for r in rows]; ob = purity(kc, yy); nl = []
for _ in range(1000): s = list(yy); random.shuffle(s); nl.append(purity(kc, s))
print("purity on mode-preserving plant: %.3f null %.3f p=%.3f (blind expected)" % (ob, sum(nl)/1000, (sum(n >= ob for n in nl)+1)/1001))
