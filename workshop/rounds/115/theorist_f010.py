"""Round 115 (theorist): HH^0..2 on the n = 9 F-010 pair 3060000 vs 3304000 (equal Coxeter key, inequivalent by q_S / Ladkani Cor 3.15, E-057, F-047) and 2400230."""
import sys
sys.path.insert(0, '.'); sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
exec(open('workshop/rounds/115/theorist_hh.py').read().split("if __name__")[0])
for s in ('3060000', '3304000', '2400230'):
    a = nk.LinearNakayamaAlgebra(9, s)
    print(s, 'key', str(search._coxeterKeyOrNone(a))[:60], hh(sorted(a.quiver.nodes), sorted(a.quiver.edges(keys=True)), procedure.relationsFrom(a), top=6))
