"""T5 / round 035: Hom_K(T,T[-1]) for the AI 2.32(b) mutation T = D + N at a vertex v, computed directly (own code, mod p linear algebra,
independent of the repo's mutation/tiltingPlus).  Right modules, paths left to right, Hom(P_y,P_x) = e_xAe_y by left multiplication.
N = [D' -g-> P_v] (degrees 0,1), D' = sum_b P_t(b), g = left mult by arrows b out of v; D = sum_{j!=v} P_j.
Prints dim Hom(T,T[m]) m=-1,0,1, the split Hom(N,D[-1]) (= sum J_i), Hom(N,N[-1]), and the formula check
  Hom(N,N[-1]) = {(y_b): y_b in J_t(b), sum_b b y_b = 0}   (J_i = {x in e_iAe_v: x b = 0 all out-arrows b}).
Algebra = path algebra / <relations> / <paths of length >= L>.  usage: python scholar_hom_nn.py"""
import itertools, sys
P=32003
def inv(a): return pow(a,P-2,P)
class Alg:
    def __init__(s,verts,arrows,rels,L):
        s.verts=verts; s.arr=dict(arrows)   # name -> (tail, head)
        paths=[(('e',x),) for x in verts]  # trivial path as ('e',x)
        # represent a path as (tail, head, tuple of arrow names)
        s.paths=[(x,x,()) for x in verts]
        cur=list(s.paths)
        for _ in range(L-1):
            nxt=[(t,s.arr[a][1],p+(a,)) for (t,h,p) in cur for a,(t2,h2) in s.arr.items() if t2==h]
            s.paths+=nxt; cur=nxt
        s.idx={p:i for i,p in enumerate(s.paths)}
        gens=[]  # ideal spanning vectors
        for r in rels:  # r: dict path-tuple -> coeff ; all terms must have same tail/head
            for (t,h,p) in s.paths:
                pass
        def pth(p):  # arrow-name tuple -> path triple
            t=s.arr[p[0]][0]; h=s.arr[p[-1]][1]; return (t,h,tuple(p))
        s.pth=pth
        for r in rels:
            terms=[(pth(p),c) for p,c in r.items()]
            t0,h0=terms[0][0][0],terms[0][0][1]
            for (u,uh,up) in s.paths:
                for (w,wh,wp) in s.paths:
                    if uh==t0 and h0==w:
                        v={}
                        for (pp,c) in terms:
                            q=(u,wh,up+pp[2]+wp)
                            if q in s.idx: v[s.idx[q]]=(v.get(s.idx[q],0)+c)%P
                        if v: gens.append(v)
        # echelon form of ideal
        s.piv={}  # pivot col -> row dict (normalised)
        for g in gens: s._add(g)
        s.basis=[i for i in range(len(s.paths)) if i not in s.piv]
        s.bpos={i:k for k,i in enumerate(s.basis)}
    def _red(s,v):
        v=dict(v)
        for c in sorted(list(v)):
            pass
        changed=True
        while True:
            cs=[c for c in v if c in s.piv and v[c]%P]
            if not cs: break
            c=cs[0]; f=v[c]
            for k,x in s.piv[c].items(): v[k]=(v.get(k,0)-f*x)%P
            v={k:x for k,x in v.items() if x%P}
        return {k:x%P for k,x in v.items() if x%P}
    def _add(s,g):
        g=s._red(g)
        if not g: return
        c=min(g); f=inv(g[c]); g={k:x*f%P for k,x in g.items()}
        for k in list(s.piv):  # keep reduced
            if c in s.piv[k]:
                f2=s.piv[k][c]; r=s.piv[k]
                for kk,x in g.items(): r[kk]=(r.get(kk,0)-f2*x)%P
                s.piv[k]={a:b for a,b in r.items() if b%P}
        s.piv[c]=g
    def mul(s,x,y):  # x,y: dict basis-path-index -> coeff (elements); product x then y
        out={}
        for i,a in x.items():
            ti,hi,pi=s.paths[i]
            for j,b in y.items():
                tj,hj,pj=s.paths[j]
                if hi!=tj: continue
                q=(ti,hj,pi+pj)
                if q in s.idx: out[s.idx[q]]=(out.get(s.idx[q],0)+a*b)%P
        return s._red(out)
    def eAe(s,i,j):  # basis of e_i A e_j : paths from i to j
        return [k for k in s.basis if s.paths[k][0]==i and s.paths[k][1]==j]
def rank(rows,ncols):
    rows=[dict(r) for r in rows]; rk=0; used=[]
    M=[ [r.get(c,0)%P for c in range(ncols)] for r in rows]
    r0=0
    for c in range(ncols):
        pr=None
        for r in range(r0,len(M)):
            if M[r][c]: pr=r;break
        if pr is None: continue
        M[r0],M[pr]=M[pr],M[r0]; f=inv(M[r0][c]); M[r0]=[x*f%P for x in M[r0]]
        for r in range(len(M)):
            if r!=r0 and M[r][c]:
                g=M[r][c]; M[r]=[(x-g*y)%P for x,y in zip(M[r],M[r0])]
        r0+=1
    return r0
class Cx:  # complex of projectives: terms[k] = list of vertices; d[k] = dict (a,b)-> element for X^k summand a -> X^{k+1} summand b ; map P_a->P_b is x in e_bAe_a
    def __init__(s,terms,d): s.terms=terms; s.d=d
def homdim(A,X,Y,m):
    """dim H^m of Hom complex Hom^*(X,Y): Hom^m = prod_k Hom(X^k,Y^{k+m}); (Df)=d_Y f - (-1)^m f d_X."""
    def coords(m):
        co=[]
        for k in X.terms:
            if k+m in Y.terms:
                for ia,a in enumerate(X.terms[k]):
                    for ib,b in enumerate(Y.terms[k+m]):
                        for e in A.eAe(b,a): co.append((k,ia,ib,e))
        return co
    def D(m):
        src=coords(m); tgt=coords(m+1); ti={c:i for i,c in enumerate(tgt)}
        rows=[]
        for (k,ia,ib,e) in src:
            out={}
            f={e:1}  # f: X^k[ia] -> Y^{k+m}[ib], element e in e_bAe_a
            # d_Y f : X^k -> Y^{k+m+1}
            if k+m in Y.d:
                for (b,c),x in Y.d[k+m].items():
                    if b==ib:  # composite P_a -> P_b -> P_c : x f
                        for kk,v in A.mul(x,f).items(): out[(k,ia,c,kk)]=(out.get((k,ia,c,kk),0)+v)%P
            # f d_X : X^{k-1} -> Y^{k+m}
            if k-1 in X.d:
                for (a0,a),x in X.d[k-1].items():
                    if a==ia:
                        for kk,v in A.mul(f,x).items(): out[(k-1,a0,ib,kk)]=(out.get((k-1,a0,ib,kk),0)-(1 if m%2==0 else -1)*v)%P
            rows.append({ti[c]:v for c,v in out.items() if v%P})
        return rows,len(src),len(tgt)
    r,ns,nt=D(m); rk=rank(r,nt) if nt else 0
    r2,ns2,nt2=D(m-1); rk2=rank(r2,nt2) if nt2 else 0
    return ns-rk-rk2
def mutation_complexes(A,v):
    outs=[(b,A.arr[b][1]) for b in A.arr if A.arr[b][0]==v]
    def el(b): return {A.idx[A.pth((b,))]:1}
    N=Cx({0:[t for b,t in outs],1:[v]},{0:{(i,0):el(b) for i,(b,t) in enumerate(outs)}})
    Dv=[j for j in A.verts if j!=v]
    D=Cx({0:Dv},{})
    # T = D + N: concatenate summands
    T=Cx({0:Dv+[t for b,t in outs],1:[v]},{0:{(len(Dv)+i,0):el(b) for i,(b,t) in enumerate(outs)}})
    return N,D,T,outs
def report(name,A,v):
    N,D,T,outs=mutation_complexes(A,v)
    h=lambda X,Y,m: homdim(A,X,Y,m)
    # J_i
    J={}
    for i in A.verts:
        if i==v: continue
        base=A.eAe(i,v); rows=[]
        for x in base:
            r={}
            for bi,(b,t) in enumerate(outs):
                for k,val in A.mul({x:1},{A.idx[A.pth((b,))]:1}).items(): r[(bi,k)]=val
            rows.append(r)
        cols=sorted({c for r in rows for c in r}); ci={c:i for i,c in enumerate(cols)}
        J[i]=len(base)-(rank([{ci[c]:val for c,val in r.items()} for r in rows],len(cols)) if cols else 0)
    sJ=sum(J.values())
    # predicted Hom(N,N[-1]) = {(y_b): y_b in J_t(b), sum_b b y_b = 0}
    cand=[]  # basis of J_t(b) as vectors
    rowsN=[];
    pred=None
    # direct: kernel of x->(y_b) ... build full space = sum_b e_t(b)Ae_v with conditions y_b b'=0, sum b y_b=0
    coords=[(bi,e) for bi,(b,t) in enumerate(outs) for e in A.eAe(t,v)]
    conds=[]
    rws=[]
    for (bi,e) in coords:
        r={}
        for bj,(b2,t2) in enumerate(outs):
            for k,val in A.mul({e:1},{A.idx[A.pth((b2,))]:1}).items(): r[('y',bi,bj,k)]=val
        for k,val in A.mul({A.idx[A.pth((outs[bi][0],))]:1},{e:1}).items(): r[('s',k)]=val
        rws.append(r)
    cols=sorted({c for r in rws for c in r},key=str); ci={c:i for i,c in enumerate(cols)}
    predN=len(coords)-(rank([{ci[c]:val for c,val in r.items()} for r in rws],len(cols)) if cols else 0)
    res={m:h(T,T,m) for m in (-2,-1,0,1,2)}
    print("== %s  v=%s outs=%s  dimA=%d"%(name,v,outs,len(A.basis)))
    print("   J_i =",J,"  sum J =",sJ)
    print("   Hom(T,T[m]) m=-2..2:",res)
    print("   Hom(N,D[-1])=%d  Hom(N,N[-1])=%d (formula %d)  Hom(D,N[-1])=%d  Hom(N,N[1])=%d  Hom(T,T[1])=%d"%(h(N,D,-1),h(N,N,-1),predN,h(D,N,-1),h(N,N,1),res[1]))
    assert h(N,N,-1)==predN and res[-1]==h(N,D,-1)+h(N,N,-1)+h(D,N,-1)+h(D,D,-1)
    assert h(N,D,-1)==sJ and res[1]==0 and res[-2]==0
    return J,h(N,N,-1),res
if __name__=="__main__":
    # 0. E-080 acyclic control: abde = acde
    A0=Alg([1,2,3,4,5],{'a':(1,2),'b':(1,3),'c':(2,4),'d':(3,4),'e':(4,5)},[{('a','c','e'):1,('b','d','e'):-1}],6)
    report("E-080 square (acyclic control)",A0,4)
    # 1. 2-cycle, rad^2 = 0: b: v->t, c: t->v, bc = cb = 0
    A1=Alg(['v','t'],{'b':('v','t'),'c':('t','v')},[],2)
    report("C2 rad^2=0 (cyclic)",A1,'v')
    # 2. 2-cycle, rad^3 = 0, no further relation: J_t = {c}?  c b != 0 here
    A2=Alg(['v','t'],{'b':('v','t'),'c':('t','v')},[],3)
    report("C2 rad^3=0",A2,'v')
    # 3. 3-cycle v->t->x->v, rad^2=0
    A3=Alg(['v','t','x'],{'b':('v','t'),'c':('t','x'),'e':('x','v')},[],2)
    report("C3 rad^2=0",A3,'v')
    # 4. cyclic, gate-style: two paths t->v, y = c g - c' d killed by b on the right and b y = 0 on the left
    arr={'b':('v','t'),'c':('t','x'),'g':('x','v'),'c2':('t','y'),'d':('y','v')}
    rels=[{('c','g','b'):1,('c2','d','b'):-1},   # y b = 0
          {('b','c','g'):1,('b','c2','d'):-1}]    # b y = 0
    A4=Alg(['v','t','x','y'],arr,rels,6)
    report("cyclic, y=cg-c'd, yb=0=by (L=6)",A4,'v')
    # 5. as 4 but only y b = 0 (b y != 0)
    A5=Alg(['v','t','x','y'],arr,rels[:1],6)
    report("cyclic, only yb=0 (L=6)",A5,'v')
