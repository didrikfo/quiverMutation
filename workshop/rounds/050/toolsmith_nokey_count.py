"""Round 050 (toolsmith), referee item 2: count no-key nodes per level from existing pickles (seconds, no expansion).
Child ball: level size (incl. no-key nodes, which are never put in `seen`) minus number of keys first seen at that depth.
Target ball: reads the cached tball pickle; prints its structure and any no-key bookkeeping.  Usage (repo root): toolsmith_nokey_count.py tag"""
import sys, pickle, collections
tag = sys.argv[1]
m = pickle.load(open('/tmp/tsm/fr_%s.meta.pkl' % tag, 'rb')); f = pickle.load(open('/tmp/tsm/fr_%s.pkl' % tag, 'rb'))
by = collections.Counter(f['seen'].values()) if isinstance(f['seen'], dict) else None
print('levels', m['levels'], 'seenN', m['seenN'], 'seen type', type(f['seen']))
if by is None: print('seen is a set; per-depth split unavailable, total no-key =', sum(m['levels']) - m['seenN'])
else:
    tot = 0
    for d, L in enumerate(m['levels']):
        nk = L - by.get(d, 0); tot += nk; print('depth', d, 'nodes', L, 'keyed', by.get(d, 0), 'no-key', nk)
    print('no-key total depth 0..6', tot)
for c in (1, 2):
    t = pickle.load(open('/tmp/tsm/tball_c%d_d5.pkl' % c, 'rb'))
    print('class', c, 'tball tuple lens', [len(x) if hasattr(x, '__len__') else x for x in t], 'depth hist', sorted(collections.Counter(t[0].values()).items()))
    print('  types', [type(x).__name__ for x in t], 'sample t[1]', str(t[1])[:200], 't[2]', str(t[2])[:200])
