"""Round 111 toolsmith, item 3: list the 89 key-kept D=0 steps of the n=7 class-1 walk checkpoint by first-expansion index (nexp) and hash the (CA, CB) pairs.
usage (repo root): .venv/bin/python workshop/rounds/111/toolsmith_c1_13_ids.py   (needs /tmp/skeptic_c1_13_ck1.pkl, i.e. the 25-minute walk)"""
import pickle, hashlib, sys
sys.path.insert(0, '.')
seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
rows = sorted([r for r in byCart.values() if r['kept'] and not r['tp']], key=lambda r: (r['nexp'], r['CA'].tobytes(), r['CB'].tobytes()))
h = hashlib.sha256()
for r in rows: h.update(r['CA'].astype('int64').tobytes() + r['CB'].astype('int64').tobytes())
print('n steps', len(rows), 'sha256 of sorted (CA,CB) int64 bytes', h.hexdigest()[:32], 'total expansions', nexp)
print('nexp of the 89:', [r['nexp'] for r in rows])
print('keys of row:', sorted(rows[0].keys()))
