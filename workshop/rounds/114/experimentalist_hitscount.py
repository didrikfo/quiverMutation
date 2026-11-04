"""Round 114 reply: count line lengths in E-263's stored hit files (round 097), per file (the files have no trailing newline).
usage: experimentalist_hitscount.py DIR   (DIR = unpacked experimentalist_scan15_hits.tar.gz)"""
import sys, glob, collections
D = sys.argv[1]
for h in (4, 5, 6, 7):
    L = collections.Counter(); bad = []
    for f in sorted(glob.glob(D + '/experimentalist_scan15_hits_%d_*.txt' % h)):
        for l in open(f):
            l = l.strip()
            if l: L[len(l)] += 1;  bad += [l] if len(l) != 13 else []
    print('h =', h, 'lines', sum(L.values()), 'by length', dict(L), 'malformed', len(bad))
    # for malformed lines, as parsed by E-263 'ends' (digit by digit): K = ends
    ks = collections.Counter()
    for l in bad:
        m = [int(c) for c in l]; R = [(a+1, a+1+x) for a, x in enumerate(m) if x]
        n = 15
        if R:
            ks[(min(s for s, e in R)-1, n-max(e for s, e in R))] += 1
    big = [k for k in ks if min(k) >= 3 or min(k) < 0 or max(k) >= 3]
    kk = collections.Counter()
    for (a, b), c in ks.items():
        for K in (a, b):
            kk['K>=3' if K >= 3 else ('K<0' if K < 0 else 'K<3')] += c
    print('   malformed lines: E-263 ends values', dict(kk), ' sample', bad[:2])
