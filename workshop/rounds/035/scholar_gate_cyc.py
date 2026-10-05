"""Does the repo gate (simple paths only) admit the cyclic cases of scholar_hom_nn.py?  Case 4: arrows b:v>t, c:t>x, g:x>v, c2:t>y, d:y>v,
relations (c g - c2 d) b = 0 and b (c g - c2 d) = 0 (plus all paths of length >= 6 zero, not seen by the gate)."""
import sys; sys.path.insert(0,'.')
import networkx as nx
from quivermutation import procedure, arrowPaths as ap
q=nx.MultiDiGraph()
E={'b':('v','t'),'c':('t','x'),'g':('x','v'),'c2':('t','y'),'d':('y','v')}
for n,(a,b) in E.items(): q.add_edge(a,b,key=n)
ar=lambda n:(E[n][0],E[n][1],n)
r1={(ar('c'),ar('g'),ar('b')):1,(ar('c2'),ar('d'),ar('b')):-1}
r2={(ar('b'),ar('c'),ar('g')):1,(ar('b'),ar('c2'),ar('d')):-1}
print("gate (yb=0,by=0):",procedure.isMutable(q,[r1,r2],'v'))
print("gate (yb=0 only):",procedure.isMutable(q,[r1],'v'))
