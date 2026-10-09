"""Row-set identity of merged-word orbits with the 444 orbit.
Usage (repo root): .venv/bin/python workshop/rounds/015/skeptic_rowset.py N
Reads the merged words from workshop/rounds/013/skeptic_orbscan_n{N}.txt (3-letter) and workshop/rounds/014/experimentalist_orbscan4_n{N}.txt (4-letter).
S = frozenset of rows of the orbit of 444 at its middle offset (closed, asserted). For every merged word, every placement: start row in S (then its closed orbit IS S
as a row set), else walk its own orbit and report size and whether it is disjoint from S. Output workshop/rounds/015/skeptic_rowset_n{N}.txt
"""
import sys
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves
n=int(sys.argv[1]); R=freeMoves.REDUCED
def words(f):
    return [l.split()[0] for l in open(f) if len(l.split())>3 and l.split()[2]=="merged"]
w3=words("workshop/rounds/013/skeptic_orbscan_n%d.txt"%n)
w4=words("workshop/rounds/014/experimentalist_orbscan4_n%d.txt"%n)
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
o444=offs("444"); rep=freeMoves.orbitReport(n,batch._rowFor(n,"444",o444[len(o444)//2]),free=R,limit=300000)
assert rep.closed; S=frozenset(rep.rows)
# row-set check: the orbit of every 444 placement equals S
out=open("workshop/rounds/015/skeptic_rowset_n%d.txt"%n,"w")
print("n",n,"|S|",len(S),file=out)
others={}
for kind,ws in (("3",w3),("4",w4)):
    for w in ws:
        os_=offs(w); inS=0; ext=[]
        for o in os_:
            st=freeMoves._startOf(batch._rowFor(n,w,o),R)
            if st in S: inS+=1
            else:
                if st not in others:
                    r=freeMoves.orbitReport(n,batch._rowFor(n,w,o),free=R,limit=300000)
                    assert r.closed; fs=frozenset(r.rows); others[st]=fs
                    for s in fs: others.setdefault(s,fs)
                fs=others[st]; ext.append((len(fs),len(fs&S)))
        tag="IN_S" if inS==len(os_) else ("PARTIAL" if inS else "OUT")
        print(kind,w,len(os_),tag,inS,sorted(set(ext)),file=out,flush=True)
print("done",file=out)
