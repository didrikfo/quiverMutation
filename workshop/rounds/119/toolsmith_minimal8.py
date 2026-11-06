"""Round 119 (toolsmith), referee item 1: for the 8 remaining failures, is any relation in the ideal of the others?
usage: toolsmith_minimal8.py (repo root)"""
import sys, random
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/118'); sys.path.insert(0, 'workshop/rounds/091')
import skeptic_conv3 as c3
import skeptic_convention as sc
from quivermutation import arrowPaths as ap
P = sc.P
sg = lambda c: c - P if c > P // 2 else c
for kind, n, seed in [('N',5,871)] + [('N',6,s) for s in (514,1629,1806,1944,2257,2794,3277)]:
    A = c3.gen(kind, n, random.Random(seed))
    q, rels = c3.repo_input(A)
    rels = [{k: sg(v) for k, v in r.items()} for r in rels]
    lens = [sorted({len(p) for p in r}) for r in rels]
    red = []
    for j, r in enumerate(rels):
        others = rels[:j] + rels[j+1:]
        red.append(bool(ap.isInIdeal(q, others, r)))
    print(seed, 'n', n, 'nrels', len(rels), 'path lengths', lens, 'redundant(each)', red, 'minimal', not any(red))
# control: seed 1614 (signed) has a redundant length-3 relation -> must report True for it
A = c3.gen('N', 5, random.Random(1614)); q, rels = c3.repo_input(A)
rels = [{k: sg(v) for k, v in r.items()} for r in rels]
print('control 1614', [bool(ap.isInIdeal(q, rels[:j]+rels[j+1:], r)) for j, r in enumerate(rels)])
