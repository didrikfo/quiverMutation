"""Round 053 (skeptic): c1 children 5 and 13 (E-158: 13 joined only by key equality with 5): compare canonicalKey, and report (None?) and Cartan equality.
Usage (repo root): skeptic_key13.py c1.pkl"""
import sys, pickle
sys.path.insert(0, '.')
from quivermutation import fingerprint, invariants
recs = [r for r in pickle.load(open(sys.argv[1], 'rb')) if r['kind'] == 'fail']
a, b = recs[5]['childObj'], recs[13]['childObj']
ka, kb = fingerprint.canonicalKey(a), fingerprint.canonicalKey(b)
print('keys equal', ka == kb, 'non-None', ka is not None, 'labelled quiver equal', sorted(a.quiver.edges) == sorted(b.quiver.edges) if hasattr(a.quiver,'edges') else '?')
print('parent depth/v', [(recs[i].get('depth'), recs[i].get('v')) for i in (5, 13)])
