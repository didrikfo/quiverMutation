"""Tabulate experimentalist_d3table_n8c0.pkl (round 038). Usage: experimentalist_d3tab.py pkl"""
import sys, pickle
from collections import Counter, defaultdict
D = pickle.load(open(sys.argv[1], 'rb')); R = D['rows']
print('expanded', D['nalg'], 'rows saved', len(R), 'of', D['allrows'], 'distinct algs', len(D['algs']))
print('(d,J) counts:', dict(sorted(Counter((r['d'], r['j']) for r in R).items())))
for name, f in [('J!=0', lambda r: r['j']), ('d>=3,J=0', lambda r: r['d'] >= 3 and not r['j'])]:
    S = [r for r in R if f(r)]
    print('==', name, len(S), 'rows;', len({r['aid'] for r in S}), 'distinct algebras;', len({(r['aid'], r['v']) for r in S}), '(alg,v)')
    print(' by (d,j,outi):', dict(sorted(Counter((r['d'], r['j'], r['outi']) for r in S).items())))
    print(' by (d,outi,outv):', dict(sorted(Counter((r['d'], r['outi'], r['outv']) for r in S).items())))
    print(' by dimA:', dict(sorted(Counter(r['dimA'] for r in S).items())))
    print(' by level:', dict(sorted(Counter(r['level'] for r in S).items())))
# separation of d=2 from d>=3 within J!=0 : J!=0 rows with d>=3
S = [r for r in R if r['j'] and r['d'] >= 3]
print('J!=0 and d>=3:', len(S))
for r in S: print('  ', {k: r[k] for k in ('exp', 'level', 'aid', 'dimA', 'v', 'i', 'd', 'j', 'outi', 'outv', 'outi_targets', 'cartrow')})
# d=2,J!=0 rows: outi distribution
T = [r for r in R if r['j'] and r['d'] == 2]
print('d=2,J!=0:', len(T), 'outi', dict(Counter(r['outi'] for r in T)), 'outv', dict(Counter(r['outv'] for r in T)), 'algs', len({r['aid'] for r in T}))
# cartan row max entry / sum
for nm, S2 in [('d=2,J!=0', T), ('d>=3,J!=0', S), ('d>=3,J=0', [r for r in R if r['d'] >= 3 and not r['j']])]:
    print(nm, 'cartan row sum', dict(sorted(Counter(sum(r['cartrow']) for r in S2).items())), 'max entry', dict(sorted(Counter(max(r['cartrow']) for r in S2).items())),
          'n zero entries', dict(sorted(Counter(r['cartrow'].count(0) for r in S2).items())), 'dimA range', (min([r['dimA'] for r in S2] or [0]), max([r['dimA'] for r in S2] or [0])))
