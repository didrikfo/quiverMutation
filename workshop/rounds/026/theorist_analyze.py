"""Round 026 (theorist): analyse collected rows. Usage: theorist_analyze.py file.json [show]"""
import sys, json
from collections import Counter
D = json.load(open(sys.argv[1])); show = len(sys.argv) > 2
def feats(r):
    v = r['v']; rels = r['rels']; outs = r['outs']
    per = []
    for b in outs:
        per.append([rel for rel in rels if any(len(q) >= 2 and q[-1] == b and q[-2] == v for q in rel)])
    return per
def kind(rel): return ('M' if len(rel) == 1 else 'S%d' % len(rel))
rows = []
for r in D['recs']:
    per = feats(r)
    if not all(per): continue
    rows.append((r, per))
print('both-carry rows', len(rows), 'rej', sum(1 for r, _ in rows if r['kerdim']))
# candidate: number of out-arrows whose carrying relation is non-monomial (commutativity), and lengths
c = Counter()
for r, per in rows:
    sig = tuple(sorted(tuple(sorted(kind(x) for x in p)) for p in per))
    c[(sig, bool(r['kerdim']))] += 1
for k, x in sorted(c.items(), key=str): print(k, x)
if show:
    for r, per in rows:
        if r['kerdim']: print('v', r['v'], 'outs', r['outs'], 'kd', r['kerdim'], [[ [ '-'.join(map(str, q)) for q in rel] for rel in p] for p in per])
