"""Round 134 response (theorist): for every par-tagged KILL printed by theorist_audit.py, test the sufficient condition
for "KILL is a change of arrows":  write the killed binomial as p - l q with p = u a w (the kept term), q = u b w,
a, b parallel.  The change of arrows a' = a - l b turns p - l q into u a' w (an automorphism of kQ).
It carries the KILL ideal to the child ideal iff no OTHER term of any relation contains the letter a (then all other
relations are fixed).  Prints, per par KILL, the verdict and whether the condition holds.
usage (repo root):  .venv/bin/python workshop/rounds/134/theorist_parcheck.py workshop/rounds/134/theorist_audit_*.out
"""
import sys, ast, re, collections
from fractions import Fraction

tab = collections.Counter()
for f in sys.argv[1:]:
    L = open(f).read().splitlines()
    for i, line in enumerate(L):
        m = re.match(r'KILL par (\S+)', line)
        if not m: continue
        rels = ast.literal_eval(L[i + 2].split('child rels ', 1)[1].replace('Fraction', ''))
        kill = ast.literal_eval(L[i + 3].split('KILL  rels ', 1)[1].replace('Fraction', ''))
        n = next(j for j in range(len(rels)) if rels[j] != kill[j])
        p = next(iter(kill[n])); q = next(t for t in rels[n] if t != p)
        d = next(j for j in range(len(p)) if p[j] != q[j]); a = p[d]
        others = [t for j, r in enumerate(rels) for t in r if not (j == n and t in (p, q))] + [q]
        cond = not any(a in t for t in others) and p.count(a) == 1
        tab[(m.group(1), 'letter a free' if cond else 'letter a reused')] += 1
        print(f.split('_')[-1], m.group(1), 'cond' if cond else 'NO-cond', 'kept', p, 'other', q)
for k, v in sorted(tab.items()): print('  ', v, k)
