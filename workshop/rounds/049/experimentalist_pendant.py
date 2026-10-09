"""Round 049 response: for each failing key-keeping edge (u,w) in a deepreplay checkpoint, is the child w unexpanded
(still in the queue: cur[pos:] or nxt), or expanded (cur[:pos]) with no admitted child? Also: out-degree of w in the edge list.
usage: experimentalist_pendant.py CKPT"""
import sys, pickle
from collections import Counter
sys.path.insert(0, '.')
S = pickle.load(open(sys.argv[1], 'rb')); E = S['edges']
cur, pos = S['cur'], S['pos']
done = {x[3] for x in cur[:pos]}; todo = {x[3] for x in cur[pos:]} | {x[3] for x in S['nxt']}
outdeg = Counter(p for p, c, f in E); indeg = Counter(c for p, c, f in E)
t = Counter()
for p, c, f in E:
    if not f: continue
    st = 'expanded' if c in done else ('unexpanded' if c in todo else 'neither')
    t[(st, 'children=%d' % outdeg[c], 'indeg=%d' % indeg[c])] += 1
for k, v in sorted(t.items()): print(k, v)
