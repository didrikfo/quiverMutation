"""Round 087 scholar, answer to referee: Serre-duality cohomology check (adapted from the skeptic scratch script ser.py).
For all k<=K, l, i: total dim H^*(Hom(P_i, nu^k P_l)) == total dim H^*(Hom(nu^(k-1) P_l, P_i)) over GF(2^31-1).
usage (repo root): .venv/bin/python workshop/rounds/087/scholar_serre_check.py RELSTRING K"""
import sys, collections
sys.path.insert(0,'workshop/rounds/087')
import scholar_rankgrowth as S
P=S.P
def rank(M):
    M=[r[:] for r in M]; r=0
    if not M: return 0
    nc=len(M[0])
    for c in range(nc):
        piv=None
        for i in range(r,len(M)):
            if M[i][c]%P: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c],P-2,P)
        M[r]=[x*inv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%P:
                f=M[i][c]; M[i]=[(x-f*y)%P for x,y in zip(M[i],M[r])]
        r+=1
    return r
def coh(A,X,i,side):
    # side 'L': Hom(P_i, X): active terms hom(i,a); side 'R': Hom(X,P_i): hom(a,i)
    ok=(lambda a: A.hom(i,a)) if side=='L' else (lambda a: A.hom(a,i))
    byd=collections.defaultdict(list)
    for u,(t,a) in X.vert.items():
        if ok(a): byd[t].append(u)
    dims={}
    ranks={}
    for t in byd:
        src=byd[t]; 
        if side=='L': tgt=byd.get(t+1,[]); 
        else: tgt=byd.get(t-1,[])  # R: d: Hom(P_b,Pi)->Hom(P_a,Pi) for a->b ; degree goes the other way
        M=[]
        if side=='L':
            for u in src:
                M.append([X.out.get(u,{}).get(v,0) for v in tgt])
        else:
            # element in Hom(X^t,P_i)=src ; image in Hom(X^{t-1}) : for w in tgt (deg t-1) with w->u
            for u in src:
                M.append([X.out.get(w,{}).get(u,0) for w in tgt])
        ranks[t]=rank(M) if tgt else 0
    res={}
    for t in byd:
        inc = ranks.get(t-1,0) if side=='L' else ranks.get(t+1,0)
        res[t]=len(byd[t])-ranks[t]-inc
    return {t:d for t,d in res.items() if d}
rs=sys.argv[1];K=int(sys.argv[2])
A,res,cur=S.run([int(c) for c in rs],K,verbose=False)
# need all levels: rerun storing
n=A.n;lv=[]
c={}
for l in range(1,n+1):
    X=S.Cx();X.vert[0]=(0,l);c[l]=X
lv.append(c)
memo={}
for k in range(1,K+1):
    nx={}
    for l,X in lv[-1].items(): nx[l]=S.nu_step(A,X,memo,400000)[0]
    lv.append(nx)
bad=0
for k in range(1,K+1):
    for l in range(1,n+1):
        for i in range(1,n+1):
            L=coh(A,lv[k][l],i,'L')           # Hom^*(P_i, nu^k P_l)
            R=coh(A,lv[k-1][l],i,'R')         # Hom^*(nu^{k-1}P_l, P_i)  ; Serre: Hom(P_i,nu Y[t]) = D Hom(Y[t],P_i)
            # Hom(Y,P_i[s]) in degree s ; L degree t = H^t ; R: H^s of Hom complex, degree s -> matching t = -s? 
            Rs=sum(R.values());Ls=sum(L.values())
            if Rs!=Ls: bad+=1; print('mismatch',k,l,i,L,R)
print(rs,'checked k<=',K,'mismatches',bad)
