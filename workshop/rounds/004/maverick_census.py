"""H-017 census: (cords, relations) of every quipu-with-relations algebra whose
Coxeter polynomial equals that of an LNA outside a quipu class, order n.
usage: python workshop/rounds/004/maverick_census.py N [minArrows]"""
import sys, collections
from quivermutation import quipuRelations as qr, coxeterTables as ct

n = int(sys.argv[1]); minA = int(sys.argv[2]) if len(sys.argv) > 2 else 2
res = qr.search(n, minArrows=minA, statuses=(ct.NOT_QUIPU, ct.UNPLACED),
                keepPerKey=None, dedupe=False)
tab = collections.defaultdict(collections.Counter)
for key, group in res['examples'].items():
    for m in group:
        cords = sum(1 for x in m['parameters'][1] if x > 0)
        tab[key][(cords, len(m['relations']))] += 1
print("walked", res['walked'], "matched", res['matched'])
for key, c in tab.items():
    print("poly", key, "LNAs", res['examples'][key][0].get('lnas'))
    for (cd, r), v in sorted(c.items()):
        print("  cords", cd, "rels", r, "count", v, "diag", r - cd)
