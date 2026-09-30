"""Summary of key-class pairing for all words over 0,3..9, <=4 letters. usage: skeptic_scan.py n [n...]"""
import sys,itertools
sys.path.insert(0,'.')
import batch
from quivermutation import coxeterTables as ct
def info(n,c):
    keys={}
    for o in range(n):
        r=batch._rowFor(n,c,o)
        if r: keys.setdefault(ct.lnaCoxeterKey(n,tuple(r)),[]).append(o)
    gs=[g for g in keys.values() if len(g)>1]
    return sorted({n-sum(g) for g in gs}),max([len(g) for g in gs],default=1)
for n in map(int,sys.argv[1:]):
    ws=[c for L in range(1,5) for w in itertools.product('03456789',repeat=L) for c in [''.join(w)]
        if c[0]!='0' and c[-1]!='0' and '00' not in c and batch._rowFor(n,c,0)]
    res={c:info(n,c) for c in ws}
    one=[c for c,(k,m) in res.items() if len(k)==1 and m==2]
    tri=[c for c,(k,m) in res.items() if m>2]
    none=[c for c,(k,m) in res.items() if not k]
    multi=[c for c,(k,m) in res.items() if len(k)>1 and m==2]
    print(f'n={n}: {len(ws)} words; one centre & pairs only: {len(one)}; class of size>=3: {len(tri)} {tri}; no equal keys: {len(none)}; several centres: {len(multi)} {multi}')
    print('  no-pair list:',' '.join(none))
