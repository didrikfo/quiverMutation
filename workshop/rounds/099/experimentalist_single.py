"""T7 controls for pair (8:4)(10:4): one bystander (t,m), k=3, margin 2; which placements lower the overlap?
Usage: experimentalist_single.py STEPS MARGIN [--jobs J]"""
import sys, time
from multiprocessing import Pool
from quivermutation import lnaMoves as lm, overlap as ov
PAIR=[(8,4),(10,4)]
def run(a):
    b,steps,margin=a
    pat=PAIR+[b]; lo=min(s for s,_ in pat); hi=max(s+k-1 for s,k in pat)
    off=max(0,7-lo); pat=[(s+off,k) for s,k in pat]; lo+=off; hi+=off
    L=hi+margin+7; rl=lm.embedPattern(L,pat,0)
    if rl is None: return (b,None)
    t0=time.time(); r=lm.localMutationSequences(L,rl,lo,hi,steps,margin); base=ov.maxOverlap(rl)
    o=[ov.maxOverlap([int(c) for c in n]) for n in r]
    return (b,"base",base,"reached",len(r),"min",min(o) if o else None,"lowered",sum(x<base for x in o),"secs",round(time.time()-t0,1))
if __name__=="__main__":
    steps,margin=int(sys.argv[1]),int(sys.argv[2]); jobs=int(sys.argv[sys.argv.index("--jobs")+1]) if "--jobs" in sys.argv else 4
    pa={a for s,k in PAIR for a in range(s,s+k)}
    cs=[(t,m) for t in list(range(4,8))+list(range(12,17)) for m in (2,3,4) if len(set(range(t,t+m))&pa)>=2 or t in(11,12)]
    cs=[(t,m) for t in list(range(4,8))+list(range(12,16)) for m in (2,3,4) if len(set(range(t,t+m))&pa)>=2]
    print("cells",len(cs),cs,flush=True)
    if "--plan" in sys.argv: sys.exit()
    with Pool(jobs) as p:
        for r in p.imap_unordered(run,[(c,steps,margin) for c in cs]): print(*r,flush=True)
