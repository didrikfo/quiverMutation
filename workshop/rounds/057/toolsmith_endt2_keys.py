"""Round 057 (toolsmith): check that the rebuilt pickle reproduces the child/parent keys printed in 049's log (determinism of toolsmith_collect.py).  CLS=2 .venv/bin/python workshop/rounds/057/toolsmith_endt2_keys.py"""
import re, hashlib, pickle
exec(open('workshop/rounds/057/toolsmith_endt2.py').read())
def tag(alg): return hashlib.md5(repr(fingerprint.canonicalKey(alg)).encode()).hexdigest()[:6]
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail'] if 'pickle' in dir() else None
txt = open('workshop/rounds/049/toolsmith_paths_logs.txt').read().split('== ')
bad = tot = 0
for sec in txt[1:]:
    name = sec.split('\n')[0].strip()
    if not name.endswith('_c%d' % CLS): continue
    for ln in sec.split('\n'):
        m = re.match(r'(paths|parents) (\d+) key (\w+)', ln)
        if m:
            r = recs[int(m.group(2))]; t = tag(r['childObj'] if m.group(1) == 'paths' else r['parentObj']); tot += 1; bad += t != m.group(3)
print('cls', CLS, 'keys checked', tot, 'mismatch', bad)
