"""n=16, 4056: are offsets 1 and 2 (equal Coxeter key) in one orbit? does each orbit hold its own mirror / the other's?"""
import sys; sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves as fm
n=16; core='4056'
R={o:fm._startOf(tuple(batch._rowFor(n,core,o)),fm.REDUCED) for o in (1,2)}
rep={o:fm.orbitReport(n,tuple(batch._rowFor(n,core,o)),free=fm.REDUCED,limit=1500000) for o in (1,2)}
for o in (1,2):
    m=fm.mirrorRow(n,R[o])
    print('offset',o,'rows',len(rep[o].rows),rep[o].stoppedBy,'holds own mirror',m in rep[o].rows,'holds other offset start',R[3-o] in rep[o].rows,'holds mirror of other start',fm.mirrorRow(n,R[3-o]) in rep[o].rows)
print('mirror of start1 == start2 ?',fm.mirrorRow(n,R[1])==R[2])
