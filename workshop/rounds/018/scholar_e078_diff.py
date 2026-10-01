"""E-078 n=5 (a=1,b=2,c=3,d=4,e=5; abde=acde) at d=4: R C R^T, Cartan(child), and the difference; also the map kernel."""
import sys; sys.path.insert(0, '.'); sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
A = pathAlgebra.PathAlgebra(); A.add_vertices_from(range(1, 6))
for x, y in [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]: A.add_arrow(x, y)
A.add_rel([[1, 2, 4, 5], [1, 3, 4, 5]])
v = 4; verts = sorted(A.vertices()); R = rplus(A, v, verts); C = cartan(A)
ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(A, v))
print('C parent\n', C); print('R C R^T\n', R.dot(C).dot(R.T)); print('C child (rewrite)\n', cartan(ch))
print('diff\n', R.dot(C).dot(R.T) - cartan(ch)); print('tiltingPlus', tiltingPlus(A.quiver, procedure.relationsFrom(A), v))
