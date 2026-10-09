"""Round 051 (toolsmith): relabeling cost of the level-6 frontier nodes of the c1 14 ball that have no canonicalKey at cap 720, and what canonicalKey costs at larger caps.
Usage (repo root): .venv/bin/python workshop/rounds/051/toolsmith_capcost.py /tmp/tsm/fr_c1_14.pkl"""
import sys, pickle, time, collections
sys.path.insert(0, '.')
from quivermutation import fingerprint
fr = pickle.load(open(sys.argv[1], 'rb'))['frontier']
cost = collections.Counter(); nok = []
for a in fr:
    c = fingerprint.relabelingCost(a.quiver)
    if c > 720: nok.append(a); cost[c] += 1
print('frontier', len(fr), 'cost>720:', len(nok), 'by cost', sorted(cost.items()))
for cap in (720, 5040, 40320, 362880):
    t = time.time(); got = 0; n = 0
    for a in nok[:12]:
        if fingerprint.relabelingCost(a.quiver) > cap: continue
        n += 1; k = fingerprint.canonicalKey(a, cap=cap); got += k is not None
    print('cap', cap, 'nodes with cost<=cap among first 12:', n, 'keys', got, '%.2fs' % (time.time() - t))
