"""How many quipus of order n have >= 2 nonpositive directions in the Euler form 2I - A (adjacency eigenvalues >= 2)?
Where the first appear, the criterion 'pos(C+C^T) <= n-2  <=>  not a quipu class' has to fail. usage: ... NMAX"""
import sys, numpy as np, networkx as nx
from quivermutation import quipuForms as qf
for n in range(4, int(sys.argv[1]) + 1):
    tot = bad = 0; ex = None
    for par in qf.allQuipusOfOrder(n):
        g = qf.graphFromQuipuParameters(*par); A = nx.to_numpy_array(g)
        big = int((np.linalg.eigvalsh(A) > 2 - 1e-9).sum()); tot += 1
        if big >= 2:
            bad += 1; ex = ex or par
    print(n, "quipus", tot, "with >=2 adjacency eigenvalues >= 2:", bad, ex, flush=True)
