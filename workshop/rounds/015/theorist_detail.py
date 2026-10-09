"""Round 015 theorist: for the n = 8 class 2 'M' parents, where does C(child) differ from R C R^T, and what is the parent there?
Usage: theorist_detail.py   (n = 8, class 2 fixed)"""
import sys, ast
exec(open('workshop/rounds/015/theorist_cartan.py').read().split("for kind, t in lines:")[0].replace("n, c = int(sys.argv[1]), int(sys.argv[2])", "n, c = 8, 2").replace("sys.argv = ['x']", "sys.argv=['x']"))
for kind, t in lines:
    if kind == 'R': d, rels, v, g, path = t
    else: d, rels, v, path = t
    alg = classes[base][path[0]]
    for w in path[1:]:
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, w))
    verts = sorted(alg.vertices()); ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, v))
    R = rplus(alg, v, verts); C = cartan(alg); D = R.dot(C).dot(R.T) - cartan(ch)
    pos = [(verts[i], verts[j]) for i in range(len(verts)) for j in range(len(verts)) if D[i, j]]
    print(kind, 'v', v, 'diff at (row,col) vertices', pos, 'arrows', sorted((a, b) for a, b, *_ in alg.quiver.edges(keys=True)))
    P = R.dot(C).dot(R.T)
    print('  predicted has negative entry:', bool((P < 0).any()), ' child-minus-predicted entries', [(verts[i], verts[j], int(cartan(ch)[i, j]), int(P[i, j])) for i in range(len(verts)) for j in range(len(verts)) if D[i, j]])
    print('  parent C row/col of v:', [int(x) for x in C[verts.index(v)]], [int(x) for x in C[:, verts.index(v)]])
