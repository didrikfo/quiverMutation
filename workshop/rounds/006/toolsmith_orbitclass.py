"""Orbit-plus-mirror classes against Coxeter-key classes over the catalogue (T3/T8, E-060).

  .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py N [--jobs J] [--max-word 4]
      runs (or resumes) `batch.py orbits N`, reads its ledger, and for every
      core whose orbits all closed compares three partitions of the offsets:
      orbits; orbit-plus-mirror (orbits joined when one holds the mirror of an
      offset of the other); key classes (`lnaCoxeterKey`).  Prints the cores
      where the partitions differ.  Key == orbit-plus-mirror is the prefilter
      claim; "key coarser" is E-060's case.
  .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py --same-orbit
      n = 16: walks the middle pairs of 4056, 46, 3355, 3445 and intersects
      the row sets (is the 20300 one orbit?).  Several minutes.
"""
import json, sys, os, subprocess
sys.path.insert(0, '.')
import batch
from quivermutation import coxeterTables as ct, freeMoves as fm


def classes(n, word, record):
    offs = record['offsets']
    orbits = record['orbits']
    parent = {o: o for o in offs}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        parent[find(a)] = find(b)
    for orb in orbits:
        for o in orb['held'][1:]:
            union(orb['held'][0], o)
    orbitPart = _part(offs, find)
    for orb in orbits:                      # join by mirror
        for o in orb['mirrors']:
            union(orb['held'][0], o)
    mirrorPart = _part(offs, find)
    keys = {}
    for o in offs:
        keys.setdefault(ct.lnaCoxeterKey(n, tuple(batch._rowFor(n, word, o))), []).append(o)
    keyPart = frozenset(frozenset(v) for v in keys.values())
    return orbitPart, mirrorPart, keyPart


def _part(offs, find):
    g = {}
    for o in offs:
        g.setdefault(find(o), set()).add(o)
    return frozenset(frozenset(v) for v in g.values())


def fmt(p):
    return "".join("{%s}" % ",".join(map(str, sorted(s))) for s in sorted(p, key=min))


def compare(n, maxWord):
    path = "logs/orbits-n{0}-w{1}a6-o1500000.jsonl".format(n, maxWord)
    recs = [json.loads(l) for l in open(path)]
    recs = {r['unit']: r for r in recs}
    tot = eq = 0; keyCoarser = []; keyFiner = []; other = []; capped = []; none = 0
    for w, r in sorted(recs.items(), key=lambda t: batch._coreSortKey(t[0])):
        res = r['result']
        if not res['orbits']:
            none += 1; continue
        if not all(o['closed'] for o in res['orbits']):
            capped.append(w); continue
        op, mp, kp = classes(n, w, res)
        tot += 1
        if mp == kp:
            eq += 1
        else:
            # refinement tests: a partition P refines Q if each block of P sits in a block of Q
            ref = lambda P, Q: all(any(b <= c for c in Q) for b in P)
            row = (w, fmt(op), fmt(mp), fmt(kp))
            if ref(mp, kp): keyCoarser.append(row)
            elif ref(kp, mp): keyFiner.append(row)
            else: other.append(row)
    print("n = %d: %d cores placed and closed, %d no placement, %d capped %s" % (n, tot, none, len(capped), capped))
    print("  orbit+mirror == key: %d;  key coarser: %d;  key finer: %d;  incomparable: %d" %
          (eq, len(keyCoarser), len(keyFiner), len(other)))
    for name, rows in (("key coarser", keyCoarser), ("key finer", keyFiner), ("incomparable", other)):
        for w, a, b, c in rows:
            print("  %-12s %-10s orbits %s | orbit+mirror %s | key %s" % (name, w, a, b, c))


def sameOrbit():
    n = 16
    want = {'4056': (1, 2), '46': (3, 4), '3355': (4, 5), '3445': (2, 3)}
    walks = {}
    for w, offs in want.items():
        for o in offs:
            row = tuple(batch._rowFor(n, w, o))
            rep = fm.orbitReport(n, row, free=fm.REDUCED, limit=1500000)
            walks[(w, o)] = (rep, fm._startOf(row, fm.REDUCED), fm._startOf(fm.mirrorRow(n, row), fm.REDUCED))
            print(w, o, len(rep.rows), rep.closed, flush=True)
    keys = list(walks)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            A, sa, ma = walks[a]; B, sb, mb = walks[b]
            inter = len(A.rows & B.rows) if hasattr(A.rows, '__and__') else len(set(A.rows) & set(B.rows))
            print("%s vs %s: shared rows %d; a-start in b %s; a-mirror in b %s" % (a, b, inter, sa in B.rows, ma in B.rows))


if __name__ == '__main__':
    argv = sys.argv[1:]
    if argv and argv[0] == '--same-orbit':
        sameOrbit(); sys.exit(0)
    n = int(argv[0]); maxWord = 4
    jobs = argv[argv.index('--jobs') + 1] if '--jobs' in argv else '4'
    rc = batch.main(["orbits", str(n), "--jobs", jobs, "--max-word", str(maxWord)])
    if rc != 0:
        sys.exit(rc)
    compare(n, maxWord)
