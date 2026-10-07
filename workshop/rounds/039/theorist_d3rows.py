"""Round 039 (theorist): tabulate theorist_d3walk_n8c0.pkl: did the child pass the key guard, by (d, J)? Usage: theorist_d3rows.py pkl"""
import sys, pickle
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); R = D['rows']
c = Counter((r['d'] if r['d'] < 3 else '>=3', r['j'] > 0, r['passed']) for r in R)
print('expanded', D['nalg'], 'rows', len(R)); print('(d, J!=0, child passes key guard): rows')
for k, v in sorted(c.items(), key=str): print('  ', k, v)
print('J != 0 and d >= 3 rows:')
for r in R:
    if r['j'] and r['d'] >= 3: print('  exp', r['exp'], 'v', r['v'], 'i', r['i'], 'd', r['d'], 'J', r['j'], 'out(i)', r['outi'], 'targets', r['outi_targets'], 'passed', r['passed'], 'dimA', r['dimA'])
# born-at-i test: J_i is "inherited" if some out-neighbour w of i also has J_w != 0 at the same (algebra, v)
from collections import defaultdict
S = defaultdict(set); nb = {}
for r in R:
    if r['j']: S[(r['exp'], r['v'])].add(r['i'])
tab = Counter(); steps = len(S)
for r in R:
    if not r['j']: continue
    inh = any(w in S[(r['exp'], r['v'])] for w in set(r['outi_targets']))
    tab[(r['d'], r['outi'], 'inherited' if inh else 'born')] += 1
print('J != 0 steps (alg, v):', steps, ' all children fail the guard:', all(r['passed'] is False for r in R if r['j']))
print('(d, out(i), born/inherited) for J != 0 rows:')
for k, c in sorted(tab.items(), key=str): print('  ', k, c)
