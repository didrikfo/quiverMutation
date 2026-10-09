"""Does the H6 ablation have power? Same placements, but walk = rules only (no free, edges, doubles)."""
import sys
from quivermutation import lnaMoves as lm, freeMoves as fm, overlap
n=int(sys.argv[1]); a=int(sys.argv[2]); b=int(sys.argv[3])
for o in range(0,n-3):
    rl=[0]*(n-2)
    if o+1>=len(rl): break
    rl[o]=a; rl[o+1]=b
    if not lm.isAdmissible(n,rl): continue
    out=[]
    for rules in (None, lm.VERIFIED_MOVES):
        w=fm.orbitReport(n,rl,rules=rules,free=False,edges=False,doubles=False,limit=20000)
        out.append((w.stoppedBy,len(w.rows),any(overlap.isAlmostSeparate(n,list(r)) for r in w.rows)))
    print(o,out,'SAME' if out[0]==out[1] else 'DIFF',flush=True)
