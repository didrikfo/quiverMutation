"""R071 scholar: run the repo's steps 1-7 on the worked example of arXiv:2112.08129 Sect. 7 (A_{2n,3}, here n = 3: line on 6 vertices,
all paths of length 3 zero) and print arrows and relations after mutation at vertex 1, then at vertex 3 (paper's numbering after relabel).
Usage: .venv/bin/python workshop/rounds/071/scholar_paper_example.py [N]  (N = 6 default)"""
import sys
from quivermutation import lnaMoves as lm, mutation, nakayama, pathAlgebra
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
r = [3] * (N - 3) + [0] * 1 if N >= 4 else None   # rows: r_j = 3 while j+3 <= N, else 0
r = [3 if j + 3 <= N else 0 for j in range(1, N - 1)]
alg = nakayama.LinearNakayamaAlgebra(N, r)
print('row', r)
def show(a, tag):
    print(tag); print(' nodes', sorted(a.quiver.nodes)); print(' edges', sorted(a.quiver.edges(keys=True)) if a.quiver.is_multigraph() else sorted(a.quiver.edges)); print(' rels', a.rels)
show(alg, 'start')
a1 = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [1])
show(a1, 'after mutation at 1')
for v in sorted(a1.quiver.nodes):
    ok = lm._quiet(mutation.mutationIsPossibleAtVertex, a1, v)
    print(' gate at', v, ok)
a2 = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(a1), [3])
show(a2, 'after mutation at 3 (paper: square moved two places right)')
