"""Round 114 (toolsmith): look at step 8 and its neighbours (n=7 class 1). Run: .venv/bin/python workshop/rounds/114/toolsmith_s8look.py"""
import sys
sys.argv = ['x', 'list']
src = open('workshop/rounds/113/toolsmith_c1_chain.py').read().split("t0 = time.time()\nnexp, Z")[0]
exec(compile(src, 'c1', 'exec'))
nexp, Z = load_steps()
for i in (6, 8):
    r = Z[i]
    for nm in ('alg', 'ch'):
        a = r[nm]
        print(i, nm, 'edges', sorted((x, y, a.quiver.number_of_edges(x, y)) for x, y in set(a.quiver.edges())))
        print('   rels', procedure.relationsFrom(a))
    print('  CA', r['CA'].tolist()); print('  CB', r['CB'].tolist())
