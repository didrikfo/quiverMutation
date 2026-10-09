"""Compare cord-member existence (maverick_predict_L3_s*.txt, n = 8, L = 3) with structural criteria on the relation lengths.
usage: maverick_criteria.py FILE...   (repository root)"""
import sys, collections
rows = []
for f in sys.argv[1:]:
    for l in open(f):
        p = l.split()
        if len(p) == 5: rows.append((int(p[0]), p[1], int(p[2])))
n = 8
def rels(seq): return [(i + 1, int(c)) for i, c in enumerate(seq) if c != "0"]
def mx(seq): return max(int(c) for c in seq)
def interior(seq):  # a relation of >= 3 arrows (path through >= 2 interior vertices)
    return mx(seq) >= 3
def nested(seq):  # a relation whose interior contains the start of another relation
    r = rels(seq)
    return any(s < s2 < s + a for s, a in r for s2, a2 in r)
def overlap(seq):  # two relations whose supports overlap in an arrow or more
    r = rels(seq)
    return any(s < s2 < s + a for s, a in r for s2, a2 in r)
crit = {"max>=3": interior}
tab = collections.defaultdict(collections.Counter)
for k, seq, c in sorted(rows):
    for name, f in crit.items(): tab[name][(f(seq), c > 0)] += 1
for name, t in tab.items(): print(name, "(pred, has cord):", dict(t))
print("total", len(rows), "with cords", sum(c > 0 for _, _, c in rows))
bad = [(k, s, c) for k, s, c in sorted(rows) if interior(s) != (c > 0)]
print("mismatches max>=3:", len(bad)); [print(b) for b in bad[:40]]
