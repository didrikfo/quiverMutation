"""Summarise skeptic_scan_n*.txt: merged rate by letter content; aaa rows; permutation-style null on 'the 4'."""
import glob,random,re
rows=[]
for f in sorted(glob.glob("workshop/rounds/010/skeptic_scan_n*.txt")):
    n=int(re.search(r"n(\d+)",f).group(1))
    for l in open(f):
        p=l.split()
        if p[0]=="done": continue
        rows.append((n,p[0],p[2]=="merged",int(p[4])))
def rate(sel,name):
    s=[r for r in rows if sel(r)]; m=sum(r[2] for r in s); print("%-34s %3d/%3d = %.2f"%(name,m,len(s),m/max(1,len(s))))
nz=lambda r:r[3]>1   # exclude size-1 (degenerate) orbits
rate(lambda r:nz(r),"all nondegenerate words")
rate(lambda r:nz(r) and "4" in r[1],"contain a 4")
rate(lambda r:nz(r) and "4" not in r[1],"no 4")
rate(lambda r:nz(r) and "4" not in r[1] and "2" not in r[1],"no 4, no 2")
rate(lambda r:nz(r) and "4" not in r[1] and "2" in r[1],"no 4, has a 2")
for c in "23456789":
    rate(lambda r,c=c:nz(r) and c in r[1] and "4" not in r[1] and "2" not in r[1] or (c in "24" and nz(r) and c in r[1] and (c=="4" or "4" not in r[1])),"letter %s (4/2 rules apart)"%c) if False else None
print("per-letter rate over words containing the letter, excluding 4 unless letter is 4:")
for c in "2356789":
    rate(lambda r,c=c:nz(r) and c in r[1] and "4" not in r[1],"  letter %s, no 4"%c)
rate(lambda r:nz(r) and "4" in r[1],"  letter 4")
print("aaa words (n, merged):")
for a in "23456789":
    print(" ",a*3,[(r[0],r[2],r[3]) for r in rows if r[1]==a*3])
# null: pick one of the letters 3..9 at random as 'special'; under a uniform-letter null with this data, how often is some letter's aaa merged at all n with >=1 n? trivial; instead:
# rank test: for each letter c in 2..9, the rate over words containing c (excl. size-1); rank of 4.
rates={}
for c in "23456789":
    s=[r for r in rows if nz(r) and c in r[1]]; rates[c]=sum(r[2] for r in s)/len(s)
print("rate over words containing c:",{c:round(v,2) for c,v in rates.items()})
