"""Round 085 (skeptic) response, item 4: tree-isomorphism output for the two names. Run from the repository root."""
import sys; sys.path.insert(0, '.')
import networkx as nx
from quivermutation import quipuForms as qf
A = qf.graphFromQuipuParameters((1,0,1),(1,4)); B = qf.graphFromQuipuParameters((1,1,2),(1,2))
print('nodes', A.number_of_nodes(), B.number_of_nodes(), 'edges', A.number_of_edges(), B.number_of_edges())
print('degree sequences', sorted(dict(A.degree()).values()), sorted(dict(B.degree()).values()))
for nm in ('canonicalTreeForm', 'adjacencySpectrumPolynomial'):
    f = getattr(qf, nm, None)
    if f is None:
        import quivermutation as q, importlib, pkgutil
        for m in pkgutil.iter_modules(q.__path__):
            mod = importlib.import_module('quivermutation.' + m.name)
            if hasattr(mod, nm): f = getattr(mod, nm); print(nm, 'found in', m.name); break
    ra, rb = f(A), f(B); print(nm, 'equal' if ra == rb else 'DIFFERENT'); print('  ', ra, '\n  ', rb)
print('networkx is_isomorphic (undirected):', nx.is_isomorphic(nx.Graph(A), nx.Graph(B)))
