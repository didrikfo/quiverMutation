import sys
sys.path.insert(0,'.')
import batch
from quivermutation import coxeterTables as ct
core=sys.argv[1]
for n in map(int,sys.argv[2:]):
    keys={}
    for o in range(n):
        r=batch._rowFor(n,core,o)
        if r: keys.setdefault(ct.lnaCoxeterKey(n,tuple(r)),[]).append(o)
    print(n,core,list(keys.values()))
