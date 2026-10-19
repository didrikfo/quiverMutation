"""Graded dims of N-normal-words (V-words, alpha^* only as last-applied letter)
vs Gamma_red (hat p, p reduced basic path of Qbar without alpha_*), per (src,tgt,degree).
Acyclic parents only. Usage: .venv/bin/python workshop/rounds/075/theorist_thetadims.py"""
from collections import Counter
def run(name, n, arrows, i):
    # arrows: (name, s, t, deg); composition right to left; path = list in order applied
    A=[a for a in arrows if a[1]==i and a[3]==0]
    chi=lambda v:1 if v==i else 0
    # Qbar arrows
    qb=[]
    for (nm,s,t,d) in arrows:
        qb.append(('bar'+nm,s,t,d-chi(s)+chi(t),'A' if (nm,s,t,d) in A else 'p'))
    for (nm,s,t,d) in A:
        qb.append(('*'+nm,t,i,0,'star')); qb.append(('_'+nm,i,t,0,'sub'))
    def paths(cur,end,deg,lst):
        yield cur,end,deg,lst
        for a in qb:
            if a[1]!=end: continue
            if a[4]=='sub': continue            # alpha_* only first applied
            if lst and lst[-1][4]=='star': continue  # alpha^* only last applied
            yield from paths(cur,a[2],deg+a[3],lst+[a])
    gam=Counter()
    for v in range(1,n+1):
        for (s,e,d,l) in paths(v,v,0,[]):
            if l and l[0][4]=='A': continue   # not reduced (starts = first applied is alpha-bar)
            if l and l[0][4]=='sub': continue # hat kills p' alpha_*
            gam[(s,e,d)]+=1
    # V letters
    V=[]
    for (nm,s,t,d) in arrows:
        if (nm,s,t,d) in A: continue
        V.append((nm,s,t,d-chi(s)+chi(t)))
        if t==i:
            for (an,as_,at,ad) in A: V.append((an+nm,s,at,d-chi(s)+chi(t)-1))
    for (nm,s,t,d) in A: V.append(('*'+nm,t,i,0))
    def words(end,deg,lst,start):
        yield start,end,deg,lst
        for a in V:
            if a[1]!=end: continue
            if lst and lst[-1][0].startswith('*'): continue
            yield from words(a[2],deg+a[3],lst+[a],start)
    N=Counter()
    for v in range(1,n+1):
        for (s,e,d,l) in words(v,0,[],v):
            if any(x[0].startswith('*') for x in l[:-1]): continue
            N[(s,e,d)]+=1
    # kM/(X) before quotient by rho would also count words with alpha^* inside; report only the normal-word count
    print(name,'i=',i,'| N dims == Gamma_red dims:',N==gam,'| total',sum(N.values()),sum(gam.values()))
    print('  by degree N:',sorted(Counter({d:0}) and Counter(k[2] for k in N.elements()).items()))
    print('  by degree G:',sorted(Counter(k[2] for k in gam.elements()).items()))
    print('  Gamma basic paths WITHOUT alpha_* (incl. non-reduced) by degree:',sorted(Counter(p[2] for v in range(1,n+1) for p in paths(v,v,0,[])).items()))
A3=[('phi',1,2,0),('psi',2,3,0),('omega',1,3,1)]
for i in (1,2,3): run('A3/rad^2',3,A3,i)
run('A3 hereditary w/ extra deg1 arrow 1->3',3,[('phi',1,2,0),('psi',2,3,0),('om',1,3,1),('ep',2,3,1)],1)
run('A3 hereditary',3,[('phi',1,2,0),('psi',2,3,0)],2)
run('Kronecker-like parallel 1->2 (deg0,deg1) +2->3',3,[('a',1,2,0),('b',1,2,1),('c',2,3,0)],1)
run('parallel a,b deg0 1->2 + 3->1',3,[('a',1,2,0),('b',1,2,0),('c',3,1,0)],1)
