"""Round 053 (skeptic): c1 children 14 and 15 (E-161: 15 joined only by key equality with 14): canonicalKey equal and non-None?
Usage (repo root): skeptic_key15.py c1.pkl"""
import sys, pickle
sys.path.insert(0, '.')
from quivermutation import fingerprint
recs = [r for r in pickle.load(open(sys.argv[1], 'rb')) if r['kind'] == 'fail']
a, b = recs[14]['childObj'], recs[15]['childObj']
ka, kb = fingerprint.canonicalKey(a), fingerprint.canonicalKey(b)
print('keys equal', ka == kb, 'non-None', ka is not None)
print('parent depth/v', [(recs[i].get('depth'), recs[i].get('v')) for i in (14, 15)])
