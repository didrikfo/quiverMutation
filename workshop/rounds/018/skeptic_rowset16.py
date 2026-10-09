"""n-length row-set test of the 444 orbit against every 4-letter word containing a 4 (letters 1..9, nondecreasing, >= MINOFF placements).
Usage (repo root): .venv/bin/python workshop/rounds/018/skeptic_rowset16.py N [MINOFF] [LIMIT] [START K M]  (LIMIT 0 = membership only)
S = closed orbit rows of 444 at its middle offset. Per word, every placement: start row in S (IN) or not (OUT: own orbit walked, size, overlap with S).
Words 3334 and 2455 are always included. Output workshop/rounds/018/skeptic_rowset16_n{N}.txt
line: word noffs tag inS orbits(size,overlap)...   tag IN_S | OUT | PARTIAL
"""
import sys, itertools
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves
n=int(sys.argv[1]); MINOFF=int(sys.argv[2]) if len(sys.argv)>2 else 4
LIMIT=int(sys.argv[3]) if len(sys.argv)>3 else 300000
START=sys.argv[4] if len(sys.argv)>4 else "0000"   # skip words < START
K=int(sys.argv[5]) if len(sys.argv)>5 else 0; M=int(sys.argv[6]) if len(sys.argv)>6 else 1   # shard K of M (by index in word list)
R=freeMoves.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
o444=offs("444"); rep=freeMoves.orbitReport(n,batch._rowFor(n,"444",o444[len(o444)//2]),free=R,limit=300000)
assert rep.closed; S=frozenset(rep.rows)
out=open("workshop/rounds/018/skeptic_rowset16_n%d%s.txt"%(n,("_fast" if LIMIT==0 else "") + ("" if M==1 and START=="0000" else "_s%d"%K)),"w")
print("n",n,"|S|",len(S),file=out,flush=True)
others={}
words=["".join(map(str,t)) for t in itertools.combinations_with_replacement(range(1,10),4) if 4 in t]
for w in ["3334","2455"]:
    if w not in words: words.append(w)
for i,w in enumerate(words):
    if w<START or i%M!=K: continue
    os_=offs(w)
    if len(os_)<MINOFF and w not in ("3334","2455"): continue
    if not os_: print(w,0,"NONE",file=out,flush=True); continue
    inS=0; ext=[]
    for o in os_:
        row=batch._rowFor(n,w,o); st=freeMoves._startOf(row,R)
        if st in S: inS+=1
        elif LIMIT==0: ext.append((0,0,0))   # LIMIT 0: membership only, no walk of outside orbits
        else:
            if st not in others:
                r=freeMoves.orbitReport(n,row,free=R,limit=LIMIT)
                fs=frozenset(r.rows)
                for s in fs: others.setdefault(s,(fs,r.closed))
                others[st]=(fs,r.closed)
            fs,cl=others[st]; ext.append((len(fs),len(fs&S),int(cl)))
    tag="IN_S" if inS==len(os_) else ("PARTIAL" if inS else "OUT")
    print(w,len(os_),tag,inS,sorted(set(ext)),file=out,flush=True)
print("done",file=out)
