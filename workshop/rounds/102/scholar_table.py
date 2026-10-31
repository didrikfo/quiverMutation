"""Round 102 scholar: one explicit mixed-digit pair per n = 9..12 from Fosse 2305.06642 Cor 4.6 (orbit member) and Cor 4.3 (added length-2 relation), with Kupisch series.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/102/scholar_table.py"""
from quivermutation.nakayama import LinearNakayamaAlgebra as L
from quivermutation import invariants
for n in range(9, 13):
    shown = 0
    for a in L.allOfLength(n):
        R = a.relations()
        if len(R) < 2 or any(l < 3 for _, l in R) or len({l for _, l in R}) < 2 or not a.hasAlmostSeparateRelations(): continue
        orb = sorted(a.classPreservingOrbit(), key=lambda t: t.relLengths)
        if len(orb) < 4: continue
        b = [x for x in orb if x != a and len({l for _, l in x.relations()}) >= 2][0]
        same = invariants.coxeterKey(a) == invariants.coxeterKey(b)
        print("n=%d Cor4.6  %s ~ %s  Kupisch %s ~ %s  coxeterKey equal=%s  orbit size %d" % (n, a.className(), b.className(), "".join(map(str, a.kupischSeries())) if n < 10 else a.kupischSeries(), "".join(map(str, b.kupischSeries())) if n < 10 else b.kupischSeries(), same, len(orb)))
        c = None
        for s in range(n - 2):
            if a.relLengths[s] == 0:
                d = list(a.relLengths); d[s] = 2
                try: c = L(n, d)
                except ValueError: continue
                if c.hasAlmostSeparateRelations() and c._withRelations(c.relations()) is not None: break
                c = None
        if c is not None:
            print("n=%d Cor4.3  %s ~ %s (adds length-2 at vertex %d)  Kupisch %s ~ %s  coxeterKey equal=%s" % (n, a.className(), c.className(), s + 1, a.kupischSeries(), c.kupischSeries(), invariants.coxeterKey(a) == invariants.coxeterKey(c)))
        shown += 1
        if shown == 1: break
