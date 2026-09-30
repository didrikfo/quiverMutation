"""Mirror join of orbit classes for named cores (T1/T3, round 009).
  timeout 10m .venv/bin/python workshop/rounds/009/experimentalist_mirrorjoin.py N WORD [WORD...] [--limit L]
Per core: orbit classes (held offsets, size), classes after joining orbits that hold each other's mirror, Coxeter-key classes."""
import sys, time
sys.path.insert(0, '.')
import batch
from quivermutation import coxeterTables as ct
sys.path.insert(0, 'workshop/rounds/006')
from toolsmith_orbitclass import classes, fmt

def run(n, word, limit):
    t0 = time.time()
    rec = batch.orbitCensus(n, word, limit)
    if rec is None:
        print(n, word, "no placement"); return
    closed = all(o['closed'] for o in rec['orbits'])
    op, mp, kp = classes(n, word, rec)
    sizes = [(o['held'], o['size']) for o in rec['orbits']]
    print(n, word, "offsets", rec['offsets'][0], "..", rec['offsets'][-1],
          "| orbits", fmt(op), "| sizes", " ".join("%s:%d" % (",".join(map(str, h)), s) for h, s in sizes),
          "| +mirror", fmt(mp), "| key", fmt(kp), "| closed" if closed else "| CAP", "%.0fs" % (time.time() - t0), flush=True)

if __name__ == '__main__':
    a = sys.argv[1:]
    limit = 300000
    if '--limit' in a:
        i = a.index('--limit'); limit = int(a[i + 1]); del a[i:i + 2]
    n = int(a[0])
    for w in a[1:]:
        run(n, w, limit)
