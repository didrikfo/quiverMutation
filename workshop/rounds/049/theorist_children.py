"""T10(i)/T4: describe the E-151 failing children (needs a pickle from rounds/046/skeptic_collect.py).
usage: theorist_children.py FILE.pkl"""
import sys, pickle, collections
recs = [r for r in pickle.load(open(sys.argv[1], 'rb')) if r['kind'] == 'fail']
for r in recs:
    pe = collections.Counter((a, b) for a, b, k in r['pe']); ce = collections.Counter((a, b) for a, b, k in r['ce'])
    print('depth', r['depth'], 'v', r['v'], 'J', r['J'], 'tp', r['tp'], 'parent arrows', len(r['pe']), 'child arrows', len(r['ce']),
          'child parallel', sorted(e for e, m in ce.items() if m > 1), 'child dimA', sum(map(sum, r['CB'])), 'parent dimA', sum(map(sum, r['CA'])))
    print('   parent rels', r['pr'][:150]); print('   child  rels', r['cr'][:200])
