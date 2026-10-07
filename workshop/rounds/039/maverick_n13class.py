"""n = 13: all LNAs with the Coxeter key of the lone 3 at (4,5); optionally the forward orbit of each (limit) and which rows move onto the lone 3.
usage: python workshop/rounds/039/maverick_n13class.py shard nshards   (writes hits to maverick_n13class_<shard>.txt)"""
import sys, time, itertools
from quivermutation import coxeterTables as ct, nakayama
n = 13
tgt = ct.lnaCoxeterKey(n, tuple([0]*4 + [3] + [0]*8))
sh, ns = int(sys.argv[1]), int(sys.argv[2]); mx = 0
t = time.time(); hits = []; cnt = 0
for lna in nakayama.allRelationLengths(n):
    cnt += 1
    if cnt % ns != sh: continue
    if mx and cnt > mx: break
    if ct.lnaCoxeterKey(n, tuple(lna)) == tgt: hits.append(lna)
print('rows scanned', cnt, 'hits', len(hits), '%.0fs' % (time.time() - t))
open('workshop/rounds/039/maverick_n13class_%d.txt' % sh, 'w').write('\n'.join(''.join(map(str, h)) for h in hits))
