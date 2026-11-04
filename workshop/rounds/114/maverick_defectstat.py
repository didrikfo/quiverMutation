"""T6/H-017: does the defect (rels-cords) of n=9 depth-5 quipu members depend on the key class? Which features determine it?
Statistic: purity(F) = fraction of members whose defect equals the modal defect of their F-cell (1 = defect is a function of F).
Null: defect shuffled across members (cell sizes fixed), 2000 perms. Controls: planted defect=f(F) (+ noise) via same code.
usage: .venv/bin/python workshop/rounds/114/maverick_defectstat.py   (reads workshop/rounds/105/experimentalist_defect_*.json)"""
import json, glob, collections, random
random.seed(1)
rows = []
for f in sorted(glob.glob('workshop/rounds/105/experimentalist_defect_*.json')):
    d = json.load(open(f))
    for m in d['members']:
        rows.append(dict(lna=d['lna'], key=tuple(m['key']), cords=m['cords'], len=m['len'], defect=m['rels'] - m['cords']))
FEATS = {'key': ['key'], 'key+cords': ['key', 'cords'], 'key+len': ['key', 'len'], 'cords+len': ['cords', 'len'],
         'key+cords+len': ['key', 'cords', 'len'], 'lna': ['lna'], 'lna+cords+len': ['lna', 'cords', 'len']}
def purity(cells, y):
    c = collections.defaultdict(collections.Counter)
    for k, v in zip(cells, y): c[k][v] += 1
    return sum(max(x.values()) for x in c.values()) / len(y)
def test(cells, y, B=2000):
    obs = purity(cells, y); yy = list(y); ge = 0; nulls = []
    for _ in range(B):
        random.shuffle(yy); p = purity(cells, yy); nulls.append(p); ge += p >= obs
    return obs, sum(nulls) / B, (ge + 1) / (B + 1)
def run(rs, label, yfun=None):
    y = [r['defect'] for r in rs] if yfun is None else [yfun(r) for r in rs]
    print(label)
    for name, fs in FEATS.items():
        cells = [tuple(r[f] for f in fs) for r in rs]
        o, m, p = test(cells, y, 1000)
        print("  %-15s cells=%3d purity=%.3f null=%.3f p=%.4f" % (name, len(set(cells)), o, m, p))
run(rows, "OBSERVED, all 645")
run([r for r in rows if r['cords'] >= 1], "OBSERVED, cords>=1 (584)")
# planted controls
run(rows, "PLANTED defect=f(key) (YES expected)", lambda r: 1 + (r['key'][2] != 0))
def noisy(r):
    v = 1 + (r['key'][2] != 0)
    return v if random.random() > .1 else v + 1
run(rows, "PLANTED defect=f(key) with 10% noise (YES expected)", noisy)
run(rows, "PLANTED defect iid random in {1..4} (NO expected)", lambda r: random.randint(1, 4))
