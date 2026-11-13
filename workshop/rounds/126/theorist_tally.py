"""Round 126 (theorist), response: tally of the four sweep outputs theorist_sweep_{24,34,44,53}.out (n differs, so no overlap).
Usage (repo root): .venv/bin/python workshop/rounds/126/theorist_tally.py"""
import re, ast, collections
par = inst = loop = th = both = allok = 0; leads = collections.Counter(); pl_rej = 0; where = {}; shape = collections.Counter()
for f in ['24', '34', '44', '53']:
    for line in open(f'workshop/rounds/126/theorist_sweep_{f}.out'):
        if line.startswith('=== ') and 'SKIPPED' not in line: par += 1; cur = line.split(':')[0][4:]
        m = re.match(r'  i=(\d+) (.*)', line)
        if not m: continue
        s = m.group(2); inst += 1
        lp = ast.literal_eval(re.search(r'parent_loops_at_i (\[.*?\])', s).group(1))
        tp = ast.literal_eval(re.search(r'thetaprime (\[.*?\])', s).group(1))
        ld = ast.literal_eval(re.search(r"lead (\[.*?\]) g=", s).group(1))
        ok = all(k in s for k in ['d2 0 nonhom 0', 'LF True', 'cycles True', 'C True', 'H1w True', 'PI True', 'SH True', 'ENG (0, 0, 0,'])
        allok += ok; loop += bool(lp); th += bool(tp); both += bool(lp) and bool(tp)
        if tp: pl_rej += "'YES'" not in s.split('g=pathlength')[1].split('ENG')[0]
        n = int(cur[1]); i = int(m.group(1)); xp = 'x%d' % (i % n + 1)
        for t, w in ld:
            leads[w] += 1; shape[('x\'·q' if w.split('·')[0] == xp else 'a*·p' if w.split('·')[0] == 'x%d*' % i else 'OTHER')] += 1
            where.setdefault(w, set()).add(f'{cur} i={m.group(1)}')
        if lp and tp: print('loop & Theta\':', cur, 'i =', m.group(1))
print(f'parents {par}, instances {inst}, all checks pass {allok}, loop at i {loop}, Theta\' non-empty {th}, both {both}, '
      f'g = path length rejected (all lam) {pl_rej} of {th}')
print(f'leads {sum(leads.values())}, distinct {len(leads)}:', dict(leads))
print('a*·p leads:', {w: sorted(where[w]) for w in leads if '*' in w.split('·')[0]})
print('lead shapes (x\' = arrow out of t(a) = vertex i + 1):', dict(shape))
