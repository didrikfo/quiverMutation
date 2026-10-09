"""Labelled shortest path between two placements. Usage: theorist_chain.py N W1 O1 W2 O2  (rows as in batch._rowFor)"""
import sys
sys.path.insert(0,'.'); sys.path.insert(0,'workshop/rounds/013')
import batch
from theorist_path import bfs
n=int(sys.argv[1]); w1,o1,w2,o2=sys.argv[2],int(sys.argv[3]),sys.argv[4],int(sys.argv[5])
res=bfs(n,tuple(batch._rowFor(n,w1,o1)),tuple(batch._rowFor(n,w2,o2)))
print(n,w1,o1,'->',w2,o2)
if isinstance(res,str) or res is None: print(res); sys.exit()
s,p=res; print('%d steps'%len(p)); print(' ',''.join(map(str,s)))
for r,l in p: print(' ',''.join(map(str,r)),'<-',l.split(' seq')[0])
