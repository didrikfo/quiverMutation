"""Round 113 toolsmith: K0 control. For every (kind, vertex) step of the parents of the 13 steps, det of K0 classes of T (must be +-1) and of the P2 variant
(T_k replaced by a duplicate P_j, j != k): must be 0.  Reuses toolsmith_c1_chain.py definitions via 'list' mode import."""
import sys
ARGV = ['x', 'list']
src = open('workshop/rounds/113/toolsmith_c1_chain.py').read().split("t0 = time.time()")[0]
sys.argv = ARGV; exec(compile(src, 'c', 'exec'))
nexp, Z = load_steps(); real = dup = 0; r_ok = d_zero = 0
for r in Z:
    for kind in 'FR':
        a = r['alg'] if kind == 'F' else dual(r['alg'])
        for v in sorted(a.quiver.nodes):
            if not mutation.mutationIsPossibleAtVertex(a, v): continue
            real += 1; r_ok += abs(k0det(a, v)) == 1
            for j in sorted(a.quiver.nodes):
                if j != v: dup += 1; d_zero += k0det(a, v, dup=j) == 0
print('real T: K0 det +-1 in %d of %d; duplicated-summand T: det 0 in %d of %d' % (r_ok, real, d_zero, dup))
