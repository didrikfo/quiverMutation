import re,sys,collections
for f in sys.argv[1:]:
    D=collections.Counter(); P=set(); PD=collections.defaultdict(set)
    for l in open(f):
        if l.startswith('F ('):
            m=re.match(r'F \((\d+), (\d+), (\w+), (\w+), \(([^)]*)\)',l)
            d=int(m[1]); p=tuple(m[5].split(', ')[:-1]) ; D[d]+=1; P.add(p); PD[d].add(p)
    print(f,dict(D),'edges',sum(D.values()),'parents',len(P),{d:len(v) for d,v in PD.items()})
