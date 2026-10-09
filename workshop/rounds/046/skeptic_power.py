"""Round 046 (skeptic): power check of skeptic_inv.invariants.  For each key class of LNAs (and duals) of length n, how many distinct invariant
signatures (minus charpoly) occur?  >1 means the invariants separate two algebras with equal Coxeter key (which then lie in different derived classes).
Usage: skeptic_power.py n [maxclasses]"""
import sys
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
from quivermutation import nakayama as nk, pathAlgebra, search, invariants as qi
sys.path.insert(0, 'workshop/rounds/046'); import skeptic_inv as si
n = int(ARGV[1]); classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
for ci, k in enumerate(order[:int(ARGV[2]) if len(ARGV) > 2 else 99]):
    sigs = {}
    for alg in classes[k]:
        C = [[int(x) for x in r] for r in qi.cartanMatrix(alg).tolist()]
        s = si.invariants(C, mmax=5); s.pop('charpoly'); sigs.setdefault(repr(sorted(s.items())), []).append(alg)
    print('class', ci, 'size', len(classes[k]), 'distinct signatures', len(sigs), [len(v) for v in sigs.values()], flush=True)
